import numpy as np
# from scipy import interpolate
from matplotlib import pyplot as plt

m = 1001
grid1d = np.linspace(-50, 50, m)
xm, ym = np.meshgrid(grid1d, grid1d)
x = grid1d[500:]

# we want to make two Gaussians with the same fwhmeff as
# the vk(atm) and G2(total) as in deconvVK.py.
fwhmeffvk_input = 1.7307 # 7 #0.7 #
fwhmeffg2_input = 2.7589 # 17 # 1.7589

scalef = 2*np.sqrt(2*np.log(2))
sigmavk = (fwhmeffvk_input/scalef) #0.1arcsec/pixel already included in xm, ym
sigmag2 = (fwhmeffg2_input/scalef)

vk = np.exp(-(xm*xm+ym*ym)/2/sigmavk**2)
g2 = np.exp(-(xm*xm+ym*ym)/2/sigmag2**2)

#for crosschecking
vk=vk/np.sum(vk)
fwhmeffvk = np.sqrt(1/np.sum(vk**2))*0.664*0.1 #0.1arcsec/pixel
print('cross-check: fwhmeffvk = %6.4f'%fwhmeffvk)

g2=g2/np.sum(g2)
fwhmeffg2 = np.sqrt(1/np.sum(g2**2))*0.664*0.1 #0.1arcsec/pixel
print('cross-check: fwhmeff of g2: %6.4f' % fwhmeffg2)

#now we can conv
vkfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(vk), s=vk.shape))
g2fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g2), s=g2.shape))
prodfft = g2fft*vkfft
tot = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                    s=prodfft.shape))
np.savetxt('tot.txt', tot.view(float))
tot = np.absolute(tot)
tot = tot/np.sum(tot)
fwhmefftot = np.sqrt(1/np.sum(tot**2))*0.664*0.1 #0.1arcsec/pixel
print('cross-check: fwhmeff of tot: %6.4f' % fwhmefftot)
print('cross-check: fwhmeff of tot (RSS): %6.4f' % np.sqrt(fwhmeffvk**2+fwhmeffg2**2))

vkN = vk/np.max(vk)
g2N = g2/np.max(g2)
totN = tot/np.max(tot)

#proper normalization in order to reproduce the plots

f, ax1 = plt.subplots(2,2, figsize=(8, 8))
ax1[0, 0].imshow(vkN, origin='lower')
ax1[0, 0].set_title('Gau1')
ax1[0, 1].imshow(g2N, origin='lower');
ax1[0, 1].set_title('Gau2')
ax1[1, 0].imshow(totN, origin='lower')
ax1[1, 0].set_title('Gau: tot')

ax1[1, 1].semilogy(x, vkN[500, 500:],'-b', label='Gau1')
ax1[1, 1].semilogy(x,g2N[500,500:],'-r',label='Gau2')
ax1[1, 1].semilogy(x,totN[500,500:],'-g',label='Gau tot')
ax1[1, 1].set_xlim(0, 30)
ax1[1, 1].set_ylim(1e-6, 10)
ax1[1, 1].legend(loc="upper right", fontsize=10)
        
#plt.show()
plt.savefig('output/convVKGau.png')

