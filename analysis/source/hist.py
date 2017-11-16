import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('output/correlate_temporal/run-9_fwhm_fitp.txt')
run = a[:,0]
nfield = a[:, 1]
tau = a[:, 2]
tauerr = a[:, 3]
alpha = a[:, 4]
alphaerr = a[:, 5]

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8,4))
#plt.figure()

idx = nfield>0
#idx = nfield>600
#idx = nfield>100
ax0.scatter(nfield[idx], tau[idx], s=10, c='r', edgecolors='r')
ax0.grid()
ax0.set_ylabel(r'$\tau$ (minutes)', {'fontsize': 16})
ax0.set_xlabel('Number of fields', {'fontsize': 16})
ax0.set_xlim(0, 1000)

ax1.scatter(nfield[idx], alpha[idx], s=10, c='r', edgecolors='r')
ax1.grid()
ax1.set_ylabel(r'Power-law index $\beta$', {'fontsize': 16})
ax1.set_xlabel('Number of fields', {'fontsize': 16})
ax1.set_xlim(0, 1000)


#usebands = 'ugriz'
#a=np.loadtxt('SDSSdata/%s/Stripe82_r_col3.txt'%usebands, skiprows =1);
#alpha = a[:,1]
#beta = a[:, 2]
#tau = a[:, 3]
#seeing = a[:, 4]
#
#plt.subplot(2,2,1)
#plt.hist(alpha, 20) #, normed=1, histtype='stepfilled', facecolor='g', alpha=0.75)
#plt.grid()
#plt.xlabel(r'$\alpha$') # ($\lambda$ dependence)')
#
#plt.subplot(2,2,2)
#plt.hist(beta, 20) #, normed=1, histtype='stepfilled', facecolor='g', alpha=0.75)
#plt.grid()
#plt.xlabel(r'$\beta$') # (spatial correlation slope)')
#
#plt.subplot(2,2,3)
#plt.hist(tau, 20) #, normed=1, histtype='stepfilled', facecolor='g', alpha=0.75)
#plt.grid()
#plt.xlabel(r'$\tau$ (seconds)')
#
#plt.subplot(2,2,4)
#plt.hist(seeing, 20) #, normed=1, histtype='stepfilled', facecolor='g', alpha=0.75)
#plt.grid()
#plt.xlabel(r'$\theta$ (arcsec)')
#
plt.tight_layout()
# plt.show()
plt.savefig('output/taubeta.png', dpi=500)
