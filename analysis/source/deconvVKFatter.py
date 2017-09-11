import argparse
# import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import signal

def main():
    parser = argparse.ArgumentParser(
        description='----- deconvVKFatter.py ---------')
    parser.add_argument('-regTerm', dest='regTerm', default="1e-10",
                        help='regularization term in FFT deconvolution')
    parser.add_argument('-rec',
                            help='input delivered PSF is the \
                            reconstructed version from earlier failed test',
                            action='store_true')
    args = parser.parse_args()

    regTerm = float(args.regTerm)
    x = np.linspace(0, 50, 501)

    atm = np.loadtxt('atm.txt')

    # if we want to cut the atm stamp smaller, or flatten a tail
    grid1d = np.linspace(-50, 50, 1001)
    xm, ym = np.meshgrid(grid1d, grid1d)
    rm = np.sqrt(xm**2 + ym**2)
    # atm[rm>500] = 0

    # make the atm narrow by pNarrow %
    pNarrow = 50
    grid1dp = grid1d * (1+pNarrow / 100.)
    f= interpolate.interp2d(grid1d, grid1d, atm)
    atm = f(grid1dp, grid1dp)

    atm = atm / np.sum(atm)
    atmN = atm / np.max(atm)

    if args.rec:
        deliver = np.loadtxt('deliverrec.txt')
    else:
        deliver = np.loadtxt('deliver.txt')
    # deliver[rm>40] = deliver[500, 900]
    deliver = deliver / np.sum(deliver)
    deliverN = deliver / np.max(deliver)

    # now we can deconv
    atmfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(atm), s=atm.shape))
    deliverfft = np.fft.fftshift(np.fft.fft2(
        np.fft.fftshift(deliver), s=deliver.shape))
    ratiofft = deliverfft / (atmfft + regTerm)
    inst = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ratiofft),
                                                    s=ratiofft.shape)))
    # inst[rm>40] = inst[500, 900]
    # inst[rm>25] = inst[500, 750]
    inst = signal.savgol_filter(inst, 15, 4)
    
    inst = inst / np.sum(inst)
    instN = inst / np.max(inst)

    # convolve back and check
    instfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(inst), s=inst.shape))
    prodfft = instfft * atmfft
    deliverrec = np.absolute(
        np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                     s=prodfft.shape)))
    deliverrec = deliverrec / np.sum(deliverrec)
    deliverrecN = deliverrec / np.max(deliverrec)

    m1 = max(np.argmax(deliverrecN == np.max(deliverrecN), axis=0))
    m2 = max(np.argmax(deliverrecN == np.max(deliverrecN), axis=1))

    deliverrec = np.roll(deliverrec, 500 - m1, axis=0)
    deliverrec = np.roll(deliverrec, 500 - m2, axis=1)
    deliverrecN = np.roll(deliverrecN, 500 - m1, axis=0)
    deliverrecN = np.roll(deliverrecN, 500 - m2, axis=1)

    m1 = max(np.argmax(deliverrecN == np.max(deliverrecN), axis=0))
    m2 = max(np.argmax(deliverrecN == np.max(deliverrecN), axis=1))
    print('recentered deliveredrec is centered at (%d, %d)' % (m1, m2))
    print('First 5 data points in \n deliver     &     deliverrec:\n')
    for n in range(500, 505):
        print('%10.8e    %10.8e'%(
            deliverN[500, n:n+1], deliverrecN[500, n:n+1]))
    print('\n')
    
    f, ax1 = plt.subplots(2, 2, figsize=(12, 8))
    ax1[0, 0].imshow(atm, origin='lower')
    ax1[0, 0].set_title('atmosphere')
    ax1[0, 1].imshow(deliver, origin='lower')
    ax1[0, 1].set_title('Delivered')
    ax1[1, 0].semilogy(x, atmN[500, 500:], '-b', label='atmosphere')
    ax1[1, 0].semilogy(x, deliverN[500, 500:], '-r', label='delivered PSF')
    ax1[1, 0].semilogy(x, instN[500, 500:], '-m.', label='instrument PSF')
    ax1[1, 0].semilogy(x, deliverrecN[m1, m2:m2 + 501],
                       '-.g', label='delivered reconstructed.')
    ax1[1, 0].set_xlim(0, 50)
    ax1[1, 0].set_ylim(1e-7, 10)
    ax1[1, 0].legend(loc="upper right", fontsize=10)
    ax1[1, 0].set_title('Log plot (atm is %2.0f%% narrower)'%pNarrow)
    
    ax1[1, 1].plot(x, atmN[500, 500:], '-b', label='atmosphere')
    ax1[1, 1].plot(x, deliverN[500, 500:], '-r', label='delivered PSF')
    ax1[1, 1].plot(x, instN[500, 500:], '-m.', label='instrument PSF')
    ax1[1, 1].plot(x, deliverrecN[m1, m2:m2 + 501],
                   '-.g', label='delivered reconstructed')
    ax1[1, 1].set_xlim(0, 5)
    ax1[1, 1].legend(loc="upper right", fontsize=10)
    ax1[1, 1].set_title('Linear plot (atm is %2.0f%% narrower)'%pNarrow)
    ax1[1, 1].set_ylim(0, 1)

    # plt.show()
    plt.savefig('output/deconvATM.png')


if __name__ == "__main__":
    main()


