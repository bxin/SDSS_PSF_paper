import numpy as np
# from scipy import stats
from matplotlib import pyplot as plt

# How you can ruin the test results?
#          change sigma2 to 9.5
#          increase image size to 120
#          calculate psf12 using Gaussian equation

sigma1 = 4.0
sigma2 = 10.0
m = 100

# turn on/off the "bad" behavior 
showBad = 1
# use/don't use the "regularization" fix
useRegTerm = 1
regTerm = 1E-12

if (showBad):
    m = 110   # problems become visible in the bottom right panel
    m = 120   # complete junk 

# make images of input psfs 
grid1d = np.linspace(-(m/2-1), m/2, m)
xm, ym = np.meshgrid(grid1d, grid1d)
psf1 = np.exp(-(xm*xm+ym*ym)/2/sigma1**2)
psf2 = np.exp(-(xm*xm+ym*ym)/2/sigma2**2)
# 
psf1 = psf1/np.sum(psf1)
psf2 = psf2/np.sum(psf2)

# psf12 = convolve(psf1, psf2)
psf1fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf1), s=psf1.shape))
psf2fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf2), s=psf2.shape))
prodfft = psf2fft*psf1fft
psf12 = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft), s=prodfft.shape))
sigma12 = np.sqrt(sigma1**2+sigma2**2)
# true analytic answer: psf12 = np.exp(-(xm*xm+ym*ym)/2/sigma12**2)

# psf1_deconv = deconvolve(psf12, psf2)
psf12fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf12), s=psf12.shape))
# ratiofft = psf12fft/psf2fft
ratiofft = psf12fft/(psf2fft + useRegTerm*regTerm)

psf1_deconv = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ratiofft),
                                                    s=ratiofft.shape)))

psf1 = psf1/np.sum(psf1)
psf2 = psf2/np.sum(psf2)
psf12 = np.absolute(psf12)/np.sum(np.absolute(psf12))
fwhmeffpsf1 = np.sqrt(1/np.sum(psf1**2))*0.664
fwhmeffpsf2 = np.sqrt(1/np.sum(psf2**2))*0.664
fwhmeffpsf12 = np.sqrt(1/np.sum(psf12**2))*0.664
print('cross-check: fwhmeffpsf1 = %6.4f'%fwhmeffpsf1)
print('cross-check: fwhmeffpsf2 = %6.4f'%fwhmeffpsf2)
print('cross-check: fwhmeffpsf12 = %6.4f'%fwhmeffpsf12)
print('cross-check: fwhmeffpsf12 (RSS) = %6.4f'%np.sqrt(fwhmeffpsf1**2+fwhmeffpsf2**2))

f, axes = plt.subplots(2,2)
axes[0,0].imshow(psf1)
axes[0,0].set_title('psf1')
axes[0,1].imshow(psf2)
axes[0,1].set_title('psf2')
axes[1,0].imshow(np.real(psf12))
axes[1,0].set_title('psf12')
axes[1,1].imshow(np.real(psf1_deconv))
axes[1,1].set_title('psf1_deconv')

plt.show()
# plt.savefig('deconvVKGau2.png')
