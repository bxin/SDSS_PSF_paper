import argparse
import numpy as np
#from scipy import stats
from matplotlib import pyplot as plt

# The code generates two 2D Gaussians (psf1 and psf2), convolve them to get psf12,
# then deconvolve psf2 from psf12, and compare the result of the deconvolution to psf1.
# more details in email to ZI on 8/19/16

# How you can ruin the test results?
#          change sigma2 to 9.5
#          increase image size to 120
#          calculate psf12 using Gaussian equation

def main():
    parser = argparse.ArgumentParser(
        description='----- fwhm_lambda.py ---------')
    parser.add_argument('-plotN', help='plot normalized PSFs (max=1)',
                        action='store_true')
    parser.add_argument('-regTerm', dest='regTerm', default=1e-10, type=float,
                        help='regularization term in FFT deconvolution')
    args = parser.parse_args()
    
    m = 120 #the side length of the arrays I use for the profiles
    sigma1 = 4 # the standard deviation of the first Gaussian
    sigma2 = 9.5 # the standard deviation of the second Gaussian
    grid1d = np.linspace(-(m/2-1), m/2, m)
    xm, ym = np.meshgrid(grid1d, grid1d)
    psf1 = np.exp(-(xm*xm+ym*ym)/2/sigma1**2)
    psf2 = np.exp(-(xm*xm+ym*ym)/2/sigma2**2)
    
    #psf12 = convolve(psf1, psf2)
    psf1fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf1), s=psf1.shape))
    psf2fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf2), s=psf2.shape))
    prodfft = psf2fft*psf1fft
    psf12 = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                        s=prodfft.shape))
    # sigma12 = np.sqrt(sigma1**2+sigma2**2)
    # psf12 = np.exp(-(xm*xm+ym*ym)/2/sigma12**2)
    
    #psf1_deconv = deconvolve(psf12, psf2)
    psf12fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf12), s=psf12.shape))
    ratiofft = psf12fft/(psf2fft+args.regTerm)
    psf1_deconv = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ratiofft),
                                                        s=ratiofft.shape)))
    
    psf1 = psf1/np.sum(psf1)
    psf2 = psf2/np.sum(psf2)
    psf12 = np.absolute(psf12)/np.sum(np.absolute(psf12))
    psf1_deconv = np.absolute(psf1_deconv)/np.sum(np.absolute(psf1_deconv))
    fwhmeffpsf1 = np.sqrt(1/np.sum(psf1**2))*0.664
    fwhmeffpsf2 = np.sqrt(1/np.sum(psf2**2))*0.664
    fwhmeffpsf12 = np.sqrt(1/np.sum(psf12**2))*0.664
    fwhmeffpsf1_deconv = np.sqrt(1/np.sum(psf1_deconv**2))*0.664
    print('cross-check: fwhmeffpsf1 = %6.4f'%fwhmeffpsf1)
    print('cross-check: fwhmeffpsf2 = %6.4f'%fwhmeffpsf2)
    print('cross-check: fwhmeffpsf12 = %6.4f'%fwhmeffpsf12)
    print('cross-check: fwhmeffpsf12 (RSS) = %6.4f'%np.sqrt(fwhmeffpsf1**2+fwhmeffpsf2**2))
    print('cross-check: fwhmeffpsf1_deconv = %6.4f'%fwhmeffpsf1_deconv)
    
    f, axes = plt.subplots(2,3)
    axes[0,0].imshow(psf1)
    axes[0,0].set_title('psf1')
    axes[0,1].imshow(psf2)
    axes[0,1].set_title('psf2')
    axes[0,2].imshow(np.real(psf12))
    axes[0,2].set_title('psf12')
    axes[1,0].imshow(np.real(psf1_deconv))
    axes[1,0].set_title('psf1_deconv')
    
    mid=int(m/2)
    x=grid1d[mid:]
    
    if (not args.plotN):
        axes[1,1].semilogy(x, psf1[mid,mid:],'-b', label='psf1')
        axes[1,1].semilogy(x, psf2[mid,mid:],'-r', label='psf2')
        axes[1,1].semilogy(x, psf12[mid,mid:],'-g', label='psf12')
        axes[1,1].semilogy(x, psf1_deconv[mid,mid:],'-m', label='psf1_dc')
        axes[1,1].legend(loc="upper right", fontsize=10)
        axes[1,1].set_xlim(0, 30)
        axes[1,1].set_ylim(1e-7, 1e-2)
        
        axes[1,2].semilogy(x, psf1[mid,mid:],'-bo', label='psf1')
        axes[1,2].semilogy(x, psf2[mid,mid:],'-rx', label='psf2')
        axes[1,2].semilogy(x, psf12[mid,mid:],'-g>', label='psf12')
        axes[1,2].semilogy(x, psf1_deconv[mid,mid:],'-m*', label='psf1_dc')
        axes[1,2].set_xlim(0, 10)
        axes[1,2].set_ylim(1e-4, 1e-2)
        #plt.show()
    else:
        psf1N=psf1/np.max(psf1)
        psf2N=psf2/np.max(psf2)
        psf12N=psf12/np.max(psf12)
        psf1_deconvN = psf1_deconv/np.max(psf1_deconv)
        
        axes[1,1].semilogy(x, psf1N[mid,mid:],'-b', label='psf1')
        axes[1,1].semilogy(x, psf2N[mid,mid:],'-r', label='psf2')
        axes[1,1].semilogy(x, psf12N[mid,mid:],'-g', label='psf12')
        axes[1,1].semilogy(x, psf1_deconvN[mid,mid:],'-m', label='psf1_dc')
        axes[1,1].legend(loc="upper right", fontsize=10)
        axes[1,1].set_xlim(0, 30)
        axes[1,1].set_ylim(1e-5, 1)
        
        axes[1,2].semilogy(x, psf1N[mid,mid:],'-bo', label='psf1')
        axes[1,2].semilogy(x, psf2N[mid,mid:],'-rx', label='psf2')
        axes[1,2].semilogy(x, psf12N[mid,mid:],'-g>', label='psf12')
        axes[1,2].semilogy(x, psf1_deconvN[mid,mid:],'-m*', label='psf1_dc')
        axes[1,2].set_xlim(0, 10)
        axes[1,2].set_ylim(5e-1, 1)
    
    plt.savefig('output/deconvVKGau2.png')

if __name__ == "__main__":
    main()
