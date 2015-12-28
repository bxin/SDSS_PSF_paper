import sys
import argparse

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits

from sdssinst import sdssinst
from sdssrun import sdssrun

def main():

    parser = argparse.ArgumentParser(
        description='----- fwhm_check_eff.py ---------')
    parser.add_argument('vname', choices=('fwhmeff', 'neff'),
                        help='variable we want to look at')
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
        calc = np.zeros(txtdata.shape[0])
        if args.vname == 'fwhmeff':
            myv = txtdata[:, 4]
        else:
            myv = txtdata[:, 8]
        if (1):
            icount = 0
            for camcol in range(1, sdss.nCamcol + 1):
                print('running on camcol#%d' % camcol)
                datafile = myRun.fitsdir + "photoField-%06d-%d.fits" % (run, camcol)
                hdulist = fits.open(datafile)
                hdu1 = hdulist[1].data
                # for ifield in range(myRun.nfields-530): #for plotting tests
                for ifield in range(myRun.nfields):
                    if ifield % 100 == 0:
                        print('field No. = %d/%d' % (ifield, myRun.nfields))
                    for iBand in range(0, sdss.nBand):
                        data = hdu1[ifield]
                        sigG1 = data['psf_sigma1'][iBand]
                        sigG2 = data['psf_sigma2'][iBand]
                        b = data['psf_b'][iBand]
                        pixScale = data['pixScale'][iBand]

                        grid1d = np.linspace(-30, 30, 601)#arcsec
                        x,y = np.meshgrid(grid1d, grid1d)
                        # for fits, radius must be in pixels
                        r2 = (x*x+y*y)* (1/pixScale)**2
                        psfG1 = np.exp(-0.5 * r2 / (sigG1 * sigG1))
                        psfG2 = b * np.exp(-0.5 * r2 / (sigG2 * sigG2))
                        # note division by beta! below:
                        psfG = psfG1 + psfG2
                        # normalized to 1 at r=0 by definition
                        psfModel = psfG 
                        psfModel = psfModel/np.sum(psfModel)
                        neff = 1/np.sum(psfModel**2) /(pixScale/0.1)**2
                        if args.vname == 'neff':
                            calc[icount] = neff
                        elif args.vname == 'fwhmeff':
                            # 0.4 arcsec/pixel
                            calc[icount] = 0.664*pixScale*np.sqrt(neff)
                        icount += 1
                    
        for camcol in range(1, sdss.nCamcol + 1):
                        
            for iBand in range(0, sdss.nBand):
                ax1[iBand, camcol - 1].plot(myv, calc, 'r',marker='.', linestyle='None')
                if iBand == 0 and runcount == 1:
                    ax1[iBand, camcol - 1].set_title('camcol=%d' % camcol)
                if args.vname == 'neff':
                    if camcol == 1 and runcount ==1:
                        text = 'band: %s' % sdss.band[iBand]
                        ax1[iBand, camcol -
                            1].text(30, 16, text, fontsize=15, ha='left', va='center')
                    ax1[iBand, camcol - 1].plot([0, 70], [0, 70],
                                                color='k', linestyle='-', linewidth=2)
                elif args.vname == 'fwhmeff':
                    if camcol == 1 and runcount ==1:
                        text = 'band: %s' % sdss.band[iBand]
                        ax1[iBand, camcol -
                            1].text(1.0, 0.8, text, fontsize=15, ha='left', va='center')
                    ax1[iBand, camcol - 1].plot([0.5, 2.5], [0.5, 2.5],
                                                color='k', linestyle='-', linewidth=2)
                ax1[iBand, camcol - 1].grid()
                if iBand ==sdss.nBand-1:
                    ax1[iBand, camcol - 1].set_xlabel(args.vname)
                if camcol ==1:
                    ax1[iBand, camcol - 1].set_ylabel('My calculation')
                
        # plt.suptitle('All units in arcsec ')
        plt.tight_layout()
        # plt.show()
        plt.savefig('output/run%d_check_%s.png'%(run, args.vname))
        sys.exit()
    
if __name__ == "__main__":
    main()
