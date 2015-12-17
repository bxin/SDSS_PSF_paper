import sys
import argparse
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

from sdssinst import sdssinst
from sdssrun import sdssrun
from sdsspsf import sdsspsf

def main():
    parser = argparse.ArgumentParser(
        description='----- fwhmvk_psfwidth.py ---------')
    parser.add_argument('widthOrigin', choices=('fits','calc'), help='origin of psf_width')
    args = parser.parse_args()

    objlist = np.loadtxt('data/Stripe82RunList.dat')
    
    sdss = sdssinst()
    
    runcount = 0
    
    f, ax1 = plt.subplots(sdss.nBand, sdss.nCamcol, sharex='col',
                        sharey='row', figsize=(12, 8)) #a plot is for a run
    
    for line in objlist:
    
        run = int(line[0])
        runcount += 1
        print('-- running on run# %d (seq.# %d)---------' % (run, runcount))
        myRun = sdssrun(run)
    
        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        fwhmvk = txtdata[:, 3]
        psf_width = txtdata[:, 4]
        if args.widthOrigin == 'calc':
            psf_width[:] = 0
            icount = 0
            for camcol in range(1, sdss.nCamcol + 1):
                datafile = myRun.fitsdir + "photoField-%06d-%d.fits" % (run, camcol)
                hdulist = fits.open(datafile)
                hdu1 = hdulist[1].data
                for ifield in range(myRun.nfields):
                    for iBand in range(0, sdss.nBand):
                        psf = sdsspsf(hdu1, ifield, iBand)
                        psf_width[icount]=psf.my_psf_width
                        icount += 1
                    
        for camcol in range(1, sdss.nCamcol + 1):
                        
            for iBand in range(0, sdss.nBand):
                ax1[iBand, camcol - 1].plot(psf_width, fwhmvk, 'r')
                if iBand == 0 and runcount == 1:
                    ax1[iBand, camcol - 1].set_title('camcol=%d' % camcol)
                if camcol == 1 and runcount ==1:
                    text = 'band: %s' % sdss.band[iBand]
                    ax1[iBand, camcol -
                        1].text(1.0, 0.8, text, fontsize=15, ha='left', va='center')
                ax1[iBand, camcol - 1].plot([0.5, 2.5], [0.5, 2.5],
                                            color='k', linestyle='-', linewidth=2)
                ax1[iBand, camcol - 1].grid()
                if iBand ==sdss.nBand-1:
                    ax1[iBand, camcol - 1].set_xlabel('psf_width')
                if camcol ==1:
                    ax1[iBand, camcol - 1].set_ylabel('FWHMvK')
                
        plt.suptitle('All units in arcsec ')
        plt.tight_layout()
        #plt.show()
        plt.savefig('output/fwhmvk_psfwidth.png')
        sys.exit()
    
if __name__ == "__main__":
    main()
