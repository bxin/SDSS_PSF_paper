
import sys

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

from sdsspsf import sdsspsf

objlist = np.loadtxt('data/Stripe82RunList.dat')
rootName = "SDSSdata"

nBand = 5
nCamcol = 6
# show which band
band = {}
band[0] = "u"
band[1] = "g"
band[2] = "r"
band[3] = "i"
band[4] = "z"

f, ax1 = plt.subplots(nBand, nCamcol, sharex='col', sharey='row',figsize=(12, 8))

for line in objlist:
    run = int(line[0])
    plt.suptitle('run %d'%run)
    print('------ running on run# %d ---------'%run)
    outdir = "%s/%d/" % (rootName, run)
    for camcol in range(1, nCamcol+1):
        print('running on camcol#%d'% camcol)
        datafile = outdir + "photoField-%06d-%d.fits" % (run, camcol)
        
        hdulist = fits.open(datafile)
        alldata = hdulist[1].data
        data = alldata[0]
        for iBand in range(0, nBand):
            kk = iBand*nCamcol + camcol

            psf = sdsspsf(data, iBand)
            ax1[iBand,camcol-1].plot(psf.r, psf.LpsfG, 'b')
            ax1[iBand,camcol-1].plot(psf.r, psf.LpsfW, 'g')
            ax1[iBand,camcol-1].plot(psf.r, psf.LpsfModel, 'r')
            ax1[iBand,camcol-1].errorbar(psf.OKprofRadii, psf.OKprofile,
                                         psf.OKprofileErr, fmt='ok')

            ax1[iBand,camcol-1].set_xlim(0, 30.0)
            ax1[iBand,camcol-1].set_ylim(-6, 0.5)
            if iBand==0:
                ax1[iBand,camcol-1].set_title('camcol=%d'%camcol)
            if camcol==1:
                text = 'band: %s'%band[iBand]
                ax1[iBand,camcol-1].text(10, -2, text, fontsize=15, ha='left', va='center')

    plt.show()
    #plt.savefig('output/run94psf.png')
    sys.exit()
