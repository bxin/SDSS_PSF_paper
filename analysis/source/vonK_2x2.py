import numpy as np
import matplotlib.pyplot as plt

a30 = np.loadtxt('data/vonK1.0.txt')
a10 = np.loadtxt('data/vonK1.0_10.txt')
a100 = np.loadtxt('data/vonK1.0_100.txt')
ainf = np.loadtxt('data/Kolm1.0.txt')

x0 = int((a30.shape[0]-1)/2) #center slice
x = np.arange(x0+1) * 0.1 #.1 arcsec per pixel
r30 = a30[x0, x0:]
r10 = a10[x0, x0:]
r100 = a100[x0, x0:]
rinf = ainf[x0, x0:]

r30 = r30/np.max(r30)
r10 = r10/np.max(r10)
r100 = r100/np.max(r100)
rinf = rinf/np.max(rinf)

fig, ax1 = plt.subplots(2, 2, figsize=(8,7))

ax1[0, 0].loglog(x, r10, '--', label=r'$L_0=10$m')
ax1[0, 0].loglog(x, r30, '-',  label=r'$L_0=30$m')
ax1[0, 0].loglog(x, r100, '-.', label=r'$L_0=100$m')
ax1[0, 0].loglog(x, rinf, ':', label=r'$L_0=\infty$')

ax1[0, 0].set_xlim(0, np.max(x))
ax1[0, 0].grid()
ax1[0, 0].set_xlabel('Radius (arcsec)', {'fontsize': 15})
ax1[0, 0].set_ylabel('Surface brightness', {'fontsize': 15})
leg = ax1[0, 0].legend(loc="lower left", fontsize=11)
leg.get_frame().set_alpha(0.5)
ax1[0, 0].text(0.9, 0.9, '(a)', transform=ax1[0, 0].transAxes)

ax1[1, 0].plot(x, r10/rinf, '--', label=r'$L_0=10$m')
ax1[1, 0].plot(x, r30/rinf, '-',  label=r'$L_0=30$m')
ax1[1, 0].plot(x, r100/rinf, '-.', label=r'$L_0=100$m')
ax1[1, 0].plot(x, rinf/rinf, ':', label=r'$L_0=\infty$')
#ax1[1, 0].grid()
ax1[1, 0].set_xlabel('Radius (arcsec)', {'fontsize': 15})
ax1[1, 0].set_ylabel('Normalized surface brightness', {'fontsize': 15})
leg = ax1[1, 0].legend(loc="lower right", fontsize=11)
leg.get_frame().set_alpha(0.5)
ax1[1, 0].text(0.9, 0.9, '(c)', transform=ax1[1, 0].transAxes)
ax1[1, 0].set_ylim(0.7, 1.8)

## right panel.
b03  = np.loadtxt('data/vonK0.5.txt')
b20  = np.loadtxt('data/vonK2.0.txt')
x2 = int((b20.shape[0]-1)/2) #center slice
# r10 means FWHM=1.0.
r10 = r30 #on left panel, r30 used to mean FWHM=1.0, L0=30m
r20 = b20[x2, x2:]
r03 = b03[x0, x0:]
r20 = r20/np.max(r20)
r03 = r03/np.max(r03)

x2 = np.arange(x2+1) * 0.1 #.1 arcsec per pixel #x2 used to be a number, now a vector
ax1[0, 1].loglog(x/0.5, r03, '--', label='FWHM=0.5" (stretched to 1.0")')
ax1[0, 1].loglog(x2/2.0, r20, '-', label='FWHM=2.0" (compressed to 1.0")')
ax1[0, 1].loglog(x, r10, '-.',  label='FWHM=1.0" ')

ax1[0, 1].set_xlim(0, np.max(x))
ax1[0, 1].set_ylim(1e-8,1)
ax1[0, 1].grid()
ax1[0, 1].set_xlabel('Radius (arcsec)', {'fontsize': 15})
ax1[0, 1].set_ylabel('Surface brightness', {'fontsize': 15})
leg = ax1[0, 1].legend(loc="lower left", fontsize=9)
leg.get_frame().set_alpha(0.5)
ax1[0, 1].text(0.9, 0.9, '(b)', transform=ax1[0, 1].transAxes)

idx = np.arange(0,x0+1,2)
x01 = int(x0/2+1)
ax1[1, 1].plot(x[:x01]/0.5, r03[:x01]/r10[idx], '--', label='FWHM=0.5" (stretched to 1.0")')
idx = np.arange(0,x0*2+1,2)
ax1[1, 1].plot(x, r20[idx]/r10, '-', label='FWHM=2.0" (compressed to 1.0")')
ax1[1, 1].plot(x, r10/r10,'-.',  label='FWHM=1.0" ')
#ax1[1, 1].grid()
ax1[1, 1].set_xlabel('Radius (arcsec)', {'fontsize': 15})
ax1[1, 1].set_ylabel('Normalized surface brightness', {'fontsize': 15})
leg = ax1[1, 1].legend(loc="lower left", fontsize=9)
leg.get_frame().set_alpha(0.5)
ax1[1, 1].text(0.9, 0.9, '(d)', transform=ax1[1, 1].transAxes)

## for both panels.
plt.tight_layout()
# plt.show()
plt.savefig('output/vonK_2x2.png', dpi=500)
