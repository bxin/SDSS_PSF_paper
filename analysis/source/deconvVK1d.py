import argparse
# import sys
import numpy as np
from scipy import interpolate
from scipy import signal
from matplotlib import pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description='----- fwhm_lambda.py ---------')
    parser.add_argument('-regTerm', dest='regTerm', default="1e-10", #type=float,
                        help='regularization term in FFT deconvolution')
    args = parser.parse_args()

    regTerm = float(args.regTerm)
    m = 1001
    grid1d = np.linspace(-50, 50, m)
    xm = grid1d
    x = grid1d[500:]
    
    # we take run94, field#1 as example,
    # from: draw_1run_1field_vK_2G.py, psf.scaleR = 1.4158208768975589
    vk = np.loadtxt('data/vonK1.0.txt')
    vk1d = vk[500, :]
    scaleR = 1.4158208768975589
    scaleV = 1.0350848010744576
    
    #for crosschecking with python output
    #vkN and g2N are normalized to 1*scaleV at (0,0), vk and g2 are normalized to sum_x,y = 1
    #scaleV is because the first data point measured is not at (0, 0), but its value is 1.0000000
    vk=vk/np.sum(vk)
    fwhmeffCheck = np.sqrt(1/np.sum(vk**2))*0.664*0.1 #0.1arcsec/pixel
    fwhmeffvk = fwhmeffCheck*scaleR
    print('cross-check: fwhmeffvk = %6.4f'%fwhmeffvk)
    
    x1=grid1d*scaleR
    f=interpolate.interp1d(x1, vk1d*scaleV) #, kind='cubic')
    vk1d = f(grid1d)
    # vk[np.sqrt(xm**2+ym**2)>3] = 0
    vk1d = vk1d/np.sum(vk1d)
    vk1dN = vk1d/np.max(vk1d)*scaleV
    
    # this file was generated specifically for this test, using fwhm_check_eff.py
    # -m pdb, stop at line#83
    g2N = np.loadtxt('output/run94_psf_2G_2D.txt')
    g2N = g2N*scaleV
    g2=g2N/np.sum(g2N)
    fwhmeffg2 = np.sqrt(1/np.sum(g2**2))*0.664*0.1 #0.1arcsec/pixel
    print('cross-check: fwhmeff of g2: %6.4f' % fwhmeffg2)

    g21d = g2[500,:]
    g21d = g21d/np.sum(g21d)
    g21dN = g21d/np.max(g21d)*scaleV
    
    #now we can deconv
    # The below only takes two arguments. but, psf is a scalar...???
    psf, remainder = signal.deconvolve(g21d, vk1d)
    psf = psf/np.sum(psf)
    psfN = psf/np.max(psf)*scaleV
    
    #convolve back and check
    g2rec = signal.convolve(psf, vk1d)
    g2recN = g2rec/np.max(g2rec)*scaleV

    m1=np.argmax(g2recN==np.max(g2recN))

    f, ax1 = plt.subplots(3,4, figsize=(12, 8))
    ax1[0, 0].imshow(vk, origin='lower')
    ax1[0, 0].set_title('von Karman')
    ax1[0, 1].imshow(g2, origin='lower');
    ax1[0, 1].set_title('Double Gau.')
    ax1[0, 2].semilogy(x, vk1dN[m1:],'-b', label='vonK')
    ax1[0, 2].semilogy(x,g21dN[m1:],'-r',label='G2')
    ax1[0, 2].semilogy(x,g2recN[m1:],'-g',label='G2 rec.')
    ax1[0, 2].set_xlim(0, 30)
    ax1[0, 2].set_ylim(1e-6, 10)
    ax1[0, 2].legend(loc="upper right", fontsize=10)
    ax1[0, 3].plot(x, vk1dN[m1:],'-b', label='vonK')
    ax1[0, 3].plot(x,g21dN[m1:],'-r',label='G2')
    ax1[0, 3].plot(x,g2recN[m1:],'-g',label='G2 rec.')
    ax1[0, 3].set_xlim(0, 2)
    ax1[0, 3].legend(loc="upper right", fontsize=10)
    
    # Now we check how double G as psf performs
    sigma = 0.15 #for the core
    # parameter#1,2
    for ip in range(2):
        if ip == 0:
            A1=5e-4
            sigma1 = 2
            G1 = A1*np.exp(-xm*xm/2/sigma1**2)
        else:
            A1=6e-5
            A1=np.log10(A1)
            G1=10**(A1-(A1-(-5.8))/30*np.abs(xm))
            
        g = np.exp(-xm*xm/2/sigma**2)+G1
        
        #convolution
        new = signal.convolve(g, vk1d,'same')
        newN = new/np.max(new)*scaleV
        
        # ax1[ip+1, 0].semilogy(x, psfN[500, 500:],'-ro',label='psf deconvolved')
        ax1[ip+1, 0].semilogy(x, g[500:],'-b',label='doubleG psf')
        ax1[ip+1, 0].semilogy(x, G1[500:],'-g',label='the wider Gau.')
        ax1[ip+1, 0].set_xlim(0, 30)
        ax1[ip+1, 0].set_ylim(1e-6, 10)
        ax1[ip+1, 0].grid()
        ax1[ip+1, 0].legend(loc="upper right", fontsize=10)
        # ax1[ip+1, 1].plot(x, psfN[500, 500:],'-ro',label='psf deconvolved')
        ax1[ip+1, 1].plot(x, g[500:],'-b',label='doubleG psf')
        ax1[ip+1, 1].plot(x, G1[500:],'-g',label='the wider Gau.')
        ax1[ip+1, 1].set_xlim(0, 3)
        ax1[ip+1, 1].grid()
        ax1[ip+1, 1].legend(loc="upper right", fontsize=10)
        ax1[ip+1, 2].semilogy(x, vk1dN[500:],'-b', label='vonK')
        ax1[ip+1, 2].semilogy(x,g21dN[500:],'-r',label='G2')
        ax1[ip+1, 2].semilogy(x,newN[m1:m1+501],'-g',label='G2 rec.')
        ax1[ip+1, 2].set_xlim(0, 30)
        ax1[ip+1, 2].set_ylim(1e-6, 10)
        ax1[ip+1, 2].legend(loc="upper right", fontsize=10)
    
        ax1[ip+1, 3].plot(x, vk1dN[500:],'-b', label='vonK')
        ax1[ip+1, 3].plot(x,g21dN[500:],'-r',label='G2')
        ax1[ip+1, 3].plot(x,newN[m1:m1+501],'-g',label='G2 rec.')
        ax1[ip+1, 3].set_xlim(0, 3)
        ax1[ip+1, 3].legend(loc="upper right", fontsize=10)
            
    #plt.show()
    plt.savefig('output/deconvVK%d.png' % 94)


if __name__ == "__main__":
    main()
