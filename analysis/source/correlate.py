import numpy as np
import matplotlib.pyplot as plt

usebands = 'ugriz'
a=np.loadtxt('SDSSdata/%s/Stripe82_r_col3.txt'%usebands, skiprows =1);
alpha = a[:,1]
beta = a[:, 2]
tau = a[:, 3]
seeing = a[:, 4]

#plt.figure(figsize=(10,4.5))
plt.figure(figsize=(6,4.5))

#plt.subplot(1,2,1)
#plt.scatter(beta, tau, s=10, c='r', edgecolors='r')
#plt.grid()
#plt.xlabel(r'$\beta$') # ($\lambda$ dependence)')
#plt.ylabel(r'$\tau$ (seconds)')

#plt.subplot(1,2,2)
plt.scatter(seeing, alpha, s=10, c='r', edgecolors='r')
plt.grid()
plt.xlabel(r'FWHM$_r$ (arcsec)') # (spatial correlation slope)')
plt.ylabel(r'$\alpha$')

#The files below have the von Karman prediction curves we show on the alpha vs seeing plots.
x = np.loadtxt('data/alpha_seeing_x.txt')
y = np.loadtxt('data/alpha_seeing_y.txt')
nline = x.shape[0]
for i in range(0, nline):
    L0 = x[i,0]
    xi = x[i,1:]
    yi = y[i,1:]
    plt.plot(xi,yi)
    if L0<8: #8m
        xa = 2.7 #where to annotate
        ya = yi[np.argmax(xi<xa)]
    else:
        xa = 0.6
        ya = yi[np.argmin(xi>xa+0.05)]
    if L0<1e4:
        plt.annotate('%dm'%L0,xy=(xa, ya));
    else:
         plt.annotate('Inf',xy=(xa, ya));
        
plt.xlim([0.5, 3])
plt.ylim([-0.45, -0.1])
plt.tight_layout()
# plt.show()
plt.savefig('output/correlate.png')
