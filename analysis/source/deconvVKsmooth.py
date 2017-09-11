import argparse
# import sys
import numpy as np
from matplotlib import pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description='----- deconvVKsmooth.py ---------')
    parser.add_argument('-regTerm', dest='regTerm', default="1e-10",
                        help='regularization term in FFT deconvolution')
    args = parser.parse_args()

    regTerm = float(args.regTerm)
    x = np.linspace(0, 50, 501)

    atm = np.loadtxt('atm.txt')

    # if we want to cut the atm stamp smaller
    # grid1d = np.linspace(-50, 50, 1001)
    # xm, ym = np.meshgrid(grid1d, grid1d)
    # rm = np.sqrt(xm**2 + ym**2)
    # ns = 500
    # atm[rm>ns] = 0
    
    atm = atm / np.sum(atm)
    atmN = atm / np.max(atm)

    inst = np.loadtxt('inst.txt')
    inst = inst / np.sum(inst)
    instN = inst / np.max(inst)

    # now we can forward conv
    atmfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(atm), s=atm.shape))
    instfft = np.fft.fftshift(np.fft.fft2(
        np.fft.fftshift(inst), s=inst.shape))
    prodfft = instfft * atmfft
    deliver = np.absolute(
        np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                     s=prodfft.shape)))
    deliver = deliver / np.sum(deliver)
    deliverN = deliver / np.max(deliver)
    
    # de-convolve back and check
    deliverfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(deliver), s=deliver.shape))
    ratiofft = deliverfft / (atmfft + regTerm)
    instrec = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ratiofft),
                                                    s=ratiofft.shape)))
    instrec = instrec / np.sum(instrec)
    instrecN = instrec / np.max(instrec)

    m1 = max(np.argmax(instrecN == np.max(instrecN), axis=0))
    m2 = max(np.argmax(instrecN == np.max(instrecN), axis=1))

    f, ax1 = plt.subplots(2, 2, figsize=(12, 8))
    ax1[0, 0].imshow(atm, origin='lower')
    ax1[0, 0].set_title('atmosphere')
    ax1[0, 1].imshow(deliver, origin='lower')
    ax1[0, 1].set_title('Delivered')
    ax1[1, 0].semilogy(x, atmN[500, 500:], '-b', label='atmosphere')
    ax1[1, 0].semilogy(x, deliverN[500, 500:], '-r', label='delivered PSF')
    ax1[1, 0].semilogy(x, instN[500, 500:], '-m.', label='instrument PSF')
    ax1[1, 0].semilogy(x, instrecN[m1, m2:m2 + 501],
                       '-.g', label='instrument reconstructed.')
    ax1[1, 0].set_xlim(0, 50)
    ax1[1, 0].set_ylim(1e-6, 10)
    ax1[1, 0].legend(loc="upper right", fontsize=10)
    ax1[1, 0].set_title('Log plot')
    ax1[1, 1].plot(x, atmN[500, 500:], '-b', label='atmosphere')
    ax1[1, 1].plot(x, deliverN[500, 500:], '-r', label='delivered PSF')
    ax1[1, 1].plot(x, instN[500, 500:], '-m.', label='instrument PSF')
    ax1[1, 1].plot(x, instrecN[m1, m2:m2 + 501],
                   '-.g', label='instrument reconstructed')
    ax1[1, 1].set_xlim(0, 5)
    ax1[1, 1].legend(loc="upper right", fontsize=10)
    ax1[1, 1].set_title('Linear plot')
    
    # plt.show()
    plt.savefig('deconvSmooth.png',dpi=500)


if __name__ == "__main__":
    main()
