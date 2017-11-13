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

plt.figure(figsize=(6,4.5))

plt.semilogy(x, r10, '-.', label=r'$L_0=10$m')
plt.semilogy(x, r30, '-',  label=r'$L_0=30$m')
plt.semilogy(x, r100, ':', label=r'$L_0=100$m')
plt.semilogy(x, rinf, '--', label=r'$L_0=\infty$')

plt.grid()
plt.xlabel('Radius (arcsec)', {'fontsize': 16})
plt.ylabel('Surface brightness', {'fontsize': 16})
leg = plt.legend(loc="upper right", fontsize=16)

plt.tight_layout()
# plt.show()
plt.savefig('output/vonK.png', dpi=500)
