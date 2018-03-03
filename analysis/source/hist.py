import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('output/correlate_temporal/run-9_fwhm_fitp.txt')
run = a[:,0]
totalT = a[:, 1]*36/60
tau = a[:, 2]
tauerr = a[:, 3]
alpha = a[:, 4]
alphaerr = a[:, 5]

fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8,4))
#plt.figure()

idx = totalT>0
#idx = totalT>600*36/60
#idx = totalT>100*36/60
ax0.scatter(totalT[idx], tau[idx], s=10, c='r', edgecolors='r')
ax0.grid()
ax0.set_ylabel(r'$\tau$ (minutes)', {'fontsize': 16})
ax0.set_xlabel('Duration of run (minutes)', {'fontsize': 16})
ax0.set_xlim(0, 600)
ax0.set_ylim(-5, 30)
ax0.plot([360, 360],[-5, 30],'b--')

ax1.scatter(totalT[idx], alpha[idx], s=10, c='r', edgecolors='r')
ax1.grid()
ax1.set_ylabel(r'Power-law index $\beta$', {'fontsize': 16})
ax1.set_xlabel('Duration of run (minutes)', {'fontsize': 16})
ax1.set_xlim(0, 600)


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
