# test forward and backward convolution using 2D Gaussians
# this was later modified and renamed deconvVKGau2.py, and sent to Zeljko as a demo as the problem with dconvolution.

import numpy as np
# from scipy import interpolate
from matplotlib import pyplot as plt
# from scipy import fftpack

regTerm = 1e-6 #0 #1e-12

m = 161 #1001  # 2001
pix = 0.1
grid1d = np.linspace(-((m-1)/2+1)*pix, (m-1)/2*pix, m-1)
xm, ym = np.meshgrid(grid1d, grid1d)
mid=int((m-1)/2)
x = grid1d[mid:]

# we want to make two Gaussians with the same fwhmeff as
# the vk(atm) and G2(total) as in deconvVK.py.
fwhmeffvk_input = 1.2 # 1.7307 # 7 #
fwhmeffg2_input = 3.2554 # 1.7589 # 17 #

scalef = 2*np.sqrt(2*np.log(2))
sigmavk = (fwhmeffvk_input/scalef) #pixel size already included in xm, ym
sigmag2 = (fwhmeffg2_input/scalef)

vk = np.exp(-(xm*xm+ym*ym)/2/sigmavk**2)
g2 = np.exp(-(xm*xm+ym*ym)/2/sigmag2**2)
#g2 = np.loadtxt('tot.txt').view(complex)
#g2 = np.real(g2)
#for crosschecking
vk=vk/np.sum(vk)
fwhmeffvk = np.sqrt(1/np.sum(vk**2))*0.664*pix 
print('cross-check: fwhmeffvk = %6.4f'%fwhmeffvk)

g2r = np.real(g2)
g2r=g2r/np.sum(g2r)
fwhmeffg2 = np.sqrt(1/np.sum(g2r**2))*0.664*pix
print('cross-check: fwhmeff of g2: %6.4f' % fwhmeffg2)

#now we can deconv
vkfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(vk), s=vk.shape))
g2fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g2), s=g2.shape))
ratiofft = g2fft/(vkfft+regTerm)
psf = np.real(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(ratiofft),
                                                    s=ratiofft.shape)))
#psf = deconvolve(g2, vk)

psf = psf/np.sum(psf)
fwhmeffpsf = np.sqrt(1/np.sum(np.real(psf)**2))*0.664*0.1 #0.1arcsec/pixel
print('cross-check: fwhmeff of psf: %6.4f' % fwhmeffpsf)
print('cross-check: fwhmeff of psf (RSS): %6.4f' % np.sqrt(fwhmeffg2**2-fwhmeffvk**2))

#convolve back and check
psffft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(psf), s=psf.shape))
prodfft = psffft*vkfft
g2rec = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                    s=prodfft.shape)))

fwhmeffg2rec = np.sqrt(1/np.sum(np.real(g2rec)**2))*0.664*0.1 #0.1arcsec/pixel
print('cross-check: fwhmeff of g2rec: %6.4f' % fwhmeffg2rec)

psf = np.real(psf)
vkN = vk/np.max(vk)
g2N = g2r/np.max(g2r)
psfN = psf/np.max(psf)
g2recN = g2rec/np.max(g2rec)

#proper normalization in order to reproduce the plots

f, ax1 = plt.subplots(2,3, figsize=(12, 8))
ax1[0, 0].imshow(vkN, origin='lower')
ax1[0, 0].set_title('Gau: Atm')
ax1[0, 1].imshow(g2N, origin='lower');
ax1[0, 1].set_title('Gau: total')
ax1[0, 2].imshow(psfN, origin='lower')
ax1[0, 2].set_title('Gau: psf')
ax1[1, 0].imshow(g2recN, origin='lower');
ax1[1, 0].set_title('Gau: total rec')

ax1[1, 1].semilogy(x, vkN[mid, mid:],'-b', label='vonK')
ax1[1, 1].semilogy(x,g2N[mid,mid:],'-r',label='G2')
ax1[1, 1].semilogy(x,g2recN[mid,mid:],'g.',label='G2 rec.')
ax1[1, 1].set_xlim(0, 30)
ax1[1, 1].set_ylim(1e-6, 10)
ax1[1, 1].legend(loc="upper right", fontsize=10)
ax1[1, 2].plot(x, vkN[mid, mid:],'-b', label='vonK')
ax1[1, 2].plot(x,g2N[mid,mid:],'-r',label='G2')
ax1[1, 2].plot(x,g2recN[mid,mid:],'g.',label='G2 rec.')
ax1[1, 2].set_xlim(0, 2)
ax1[1, 2].legend(loc="upper right", fontsize=10)

#plt.show()
plt.savefig('output/deconvVKGau.png')

