import argparse
# import sys
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description='----- fwhm_lambda.py ---------')
    parser.add_argument('-plotN', help='plot normalized PSFs (max=1)',
                        action='store_true')
    parser.add_argument('-regTerm', dest='regTerm', default="1e-10", #type=float,
                        help='regularization term in FFT deconvolution')
    args = parser.parse_args()

    regTerm = float(args.regTerm)
    m = 1001
    grid1d = np.linspace(-50, 50, m)
    xm, ym = np.meshgrid(grid1d, grid1d)
    x = grid1d[500:]
    
    # we take run94, field#1 as example,
    # from: draw_1run_1field_vK_2G.py, psf.scaleR = 1.4158208768975589
    vk = np.loadtxt('data/vonK1.0.txt')
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
    f=interpolate.RectBivariateSpline(x1, x1, vk*scaleV, kx=1, ky=1)
    vk = f(grid1d, grid1d)
    # vk[np.sqrt(xm**2+ym**2)>3] = 0

    # make the atm narrow by pNarrow %
    # pNarrow = 10
    # grid1dp = grid1d * (1+pNarrow / 100.)
    # f= interpolate.interp2d(grid1d, grid1d, vk)
    # vk = f(grid1dp, grid1dp)
    
    vk = vk/np.sum(vk)
    vkN = vk/np.max(vk)*scaleV
    
    # this file was generated specifically for this test, using fwhm_check_eff.py
    # -m pdb, stop at line#83
    g2N = np.loadtxt('output/run94_psf_2G_2D.txt')
    g2N = g2N*scaleV
    g2=g2N/np.sum(g2N)
    fwhmeffg2 = np.sqrt(1/np.sum(g2**2))*0.664*0.1 #0.1arcsec/pixel
    print('cross-check: fwhmeff of g2: %6.4f' % fwhmeffg2)
    
    #now we can deconv
    vkfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(vk), s=vk.shape))
    g2fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g2), s=g2.shape))
    ratiofft = g2fft/(vkfft+regTerm)
    psf = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ratiofft),
                                                        s=ratiofft.shape)))
    # psf[np.sqrt(xm**2+ym**2)>0.4] = 0
    psf = psf/np.sum(psf)
    psfN = psf/np.max(psf)*scaleV
    
    #convolve back and check
    psffft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf), s=psf.shape))
    prodfft = psffft*vkfft
    g2rec = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                        s=prodfft.shape)))
    g2recN = g2rec/np.max(g2rec)*scaleV

    m1=max(np.argmax(g2recN==np.max(g2recN), axis=0))
    m2=max(np.argmax(g2recN==np.max(g2recN), axis=1))

    vdata = np.loadtxt('output/run94_data_points.txt')
    deliverR = vdata[0]
    deliverY = vdata[1]
    deliverE = vdata[2]
    
    # 
    # f, ax1 = plt.subplots(1,2, figsize=(12, 8))
    # ax1[0].semilogy(x, vkN[500, 500:],'-b', label='vonK')
    # ax1[0].semilogy(x,g2N[500,500:],'-r',label='G2')
    # ax1[0].semilogy(x,psfN[500,500:],'-r',label='psf')    
    # ax1[0].semilogy(x,g2recN[m1,m2:m2+501],'-g',label='G2 rec.')
    # 
    # ax1[1].plot(x, vkN[500, 500:],'-b', label='vonK')
    # ax1[1].plot(x,g2N[500,500:],'-r',label='G2')
    # ax1[1].plot(x,psfN[500,500:],'-r',label='psf')
    # ax1[1].plot(x,g2recN[m1,m2:m2+501],'-g',label='G2 rec.')
    # ax1[1].set_xlim(0, 2)
    # ax1[1].legend(loc="upper right", fontsize=10)
    # 
    # # plt.show()
    # plt.savefig('output/deconvVK%d.png' % 94)
    # sys.exit()
    #proper normalization in order to reproduce the plots
    
    f, ax1 = plt.subplots(3,4, figsize=(12, 8))
    ax1[0, 0].imshow(vk, origin='lower')
    ax1[0, 0].set_title('von Karman')
    ax1[0, 1].imshow(g2, origin='lower');
    ax1[0, 1].set_title('Double Gau.')
    ax1[0, 2].semilogy(x, vkN[500, 500:],'-b', label='vonK')
    ax1[0, 2].semilogy(x,g2N[500,500:],'-r',label='G2')
    ax1[0, 2].semilogy(x,g2recN[m1,m2:m2+501],'-g',label='G2 rec.')
    ax1[0, 2].set_xlim(0, 30)
    ax1[0, 2].set_ylim(1e-6, 10)
    ax1[0, 2].legend(loc="upper right", fontsize=10)
    ax1[0, 2].errorbar(deliverR, deliverY, deliverE, fmt='ko')
    ax1[0, 3].plot(x, vkN[500, 500:],'-b', label='vonK')
    ax1[0, 3].plot(x,g2N[500,500:],'-r',label='G2')
    ax1[0, 3].plot(x,g2recN[m1,m2:m2+501],'-g',label='G2 rec.')
    ax1[0, 3].set_xlim(0, 2)
    ax1[0, 3].legend(loc="upper right", fontsize=10)
    ax1[0, 3].errorbar(deliverR, deliverY, deliverE, fmt='ko')
    
    # Now we check how double G as psf performs
    sigma = 0.1 #for the core
    # parameter#1,2
    for ip in range(2):
        if ip == 0:
            A1=5e-4
            sigma1 = 2
            G1 = A1*np.exp(-x*x/2/sigma1**2)
            G12d = A1*np.exp(-(xm*xm+ym*ym)/2/sigma1**2)
        else:
            #A1=1.2e-5
            #sigma1 = 1200
            #G1 = A1*np.exp(-x*x/2/sigma1**2)
            #G12d = A1*np.exp(-(xm*xm+ym*ym)/2/sigma1**2)
    
            #A1=4e-5
            #G1 = A1-A1*(0.99)/20*x
            #G12d = A1-A1*(0.99)/20*np.sqrt(xm*xm+ym*ym)
            #G12d[G12d<0]=0
            
            A1=1e-5
            A1=np.log10(A1)
            x1=0
            y1=-5
            x2 = 30
            y2 = -6.9 #-6
            x3 = 15
            y3 = -6.1
            abcM = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
            abc = np.dot(np.linalg.inv(abcM), np.array([[y1] ,[y2], [y3]]))
            G1=10**(abc[0][0]*x*x+abc[1][0]*x+abc[2][0])
            rm = np.sqrt(xm*xm+ym*ym)
            G12d = 10**(abc[0][0]*rm*rm+abc[1][0]*rm+abc[2][0])
            
        g = np.exp(-x*x/2/sigma**2)+G1
        g2d = np.exp(-(xm*xm+ym*ym)/2/sigma**2)+G12d
        #m0 = int((len(grid1d)-1)/2)
        #g2d[m0, m0]= 20
        
        #convolution
        psffft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g2d), s=g2d.shape))
        prodfft = psffft*vkfft
        new = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                            s=prodfft.shape)))
        newN = new/np.max(new)*scaleV
        
        ax1[ip+1, 0].semilogy(x, psfN[500, 500:],'-ro',label='psf deconvolved')
        ax1[ip+1, 0].semilogy(x, g,'-b',label='doubleG psf')
        ax1[ip+1, 0].semilogy(x, G1,'-g',label='the wider Gau.')
        ax1[ip+1, 0].set_xlim(0, 30)
        ax1[ip+1, 0].set_ylim(1e-7, 10)
        ax1[ip+1, 0].grid()
        ax1[ip+1, 0].legend(loc="upper right", fontsize=10)
        ax1[ip+1, 1].plot(x, psfN[500, 500:],'-ro',label='psf deconvolved')
        ax1[ip+1, 1].plot(x, g,'-b',label='doubleG psf')
        ax1[ip+1, 1].plot(x, G1,'-g',label='the wider Gau.')
        ax1[ip+1, 1].set_xlim(0, 3)
        ax1[ip+1, 1].grid()
        ax1[ip+1, 1].legend(loc="upper right", fontsize=10)
        ax1[ip+1, 2].semilogy(x, vkN[500, 500:],'-b', label='vonK')
        ax1[ip+1, 2].semilogy(x,g2N[500,500:],'-r',label='G2')
        ax1[ip+1, 2].semilogy(x,newN[500,500:],'-g',label='G2 rec.')
        ax1[ip+1, 2].set_xlim(0, 30)
        ax1[ip+1, 2].set_ylim(1e-6, 10)
        ax1[ip+1, 2].legend(loc="upper right", fontsize=10)
        ax1[ip+1, 2].errorbar(deliverR, deliverY, deliverE, fmt='ko')
    
        ax1[ip+1, 3].plot(x, vkN[500, 500:],'-b', label='vonK')
        ax1[ip+1, 3].plot(x,g2N[500,500:],'-r',label='G2')
        ax1[ip+1, 3].plot(x,newN[m1,m2:m2+501],'-g',label='G2 rec.')
        ax1[ip+1, 3].set_xlim(0, 3)
        ax1[ip+1, 3].legend(loc="upper right", fontsize=10)
        ax1[ip+1, 3].errorbar(deliverR, deliverY, deliverE, fmt='ko')
        
    #plt.show()
    plt.savefig('output/deconvVK%d.png' % 94)


if __name__ == "__main__":
    main()
