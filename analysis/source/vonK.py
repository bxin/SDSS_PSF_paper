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

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8,4))

ax0.loglog(x, r10, '--', label=r'$L_0=10$m')
ax0.loglog(x, r30, '-',  label=r'$L_0=30$m')
ax0.loglog(x, r100, ':', label=r'$L_0=100$m')
ax0.loglog(x, rinf, '-.', label=r'$L_0=\infty$')

ax0.set_xlim(0, 30) #np.max(x))
ax0.grid()
ax0.set_xlabel('Radius (arcsec)', {'fontsize': 16})
ax0.set_ylabel('Surface brightness', {'fontsize': 16})
leg = ax0.legend(loc="lower left", fontsize=11)

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
ax1.loglog(x/0.5, r03, '--', label='FWHM=0.5" (stretched to 1.0")')
ax1.loglog(x, r10, '-',  label='FWHM=1.0')
ax1.loglog(x2/2.0, r20, ':', label='FWHM=2.0" (compressed to 1.0")')

ax1.set_xlim(0, 30) #np.max(x))
ax1.set_ylim(1e-8,1)
ax1.grid()
ax1.set_xlabel('Radius (arcsec)', {'fontsize': 16})
ax1.set_ylabel('Surface brightness', {'fontsize': 16})
leg = ax1.legend(loc="lower left", fontsize=9)

## for both panels.
plt.tight_layout()
# plt.show()
plt.savefig('output/vonK.png', dpi=500)
