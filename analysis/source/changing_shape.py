
import sys
import time
import argparse

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

from sdsspsf import sdsspsf


def main():
    parser = argparse.ArgumentParser(
        description='----- changing_shape.py ---------')
    parser.add_argument('yscale', choices=(
        'log', 'linear', 'logzoom'), help='yscale of the plots')
    parser.add_argument('-run', dest='irun', type=int, default=3388,
                        help='run number. default=3388, based on find_long_wild_runs.py')
    args = parser.parse_args()

    start = time.time()
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if args.irun not in objlist[:, 0]:
        print('run# %d is not in Stripe82RunList.dat\n' % args.irun)
        sys.exit()

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
    color = {}
    color[0] = "r"
    color[1] = "b"
    
    vdata = np.loadtxt('data/r_vonK_Kolm.txt',
                       unpack='True')
    radius = vdata[0]
    vonK = vdata[1]
    vonK1arcsec = np.vstack((radius, vonK))

    run = args.irun
    print('------ running on run# %d ---------' % run)
    f, ax1 = plt.subplots(nBand, nCamcol, sharex='col',
                          sharey='row', figsize=(12, 8))  # a plot is for a run

    datadir = "%s/%d/" % (rootName, run)
    fieldlist = np.zeros((nCamcol, 2))
    for camcol in range(1, nCamcol + 1):
        print('running on camcol#%d' % camcol)
        datafile = datadir + "photoField-%06d-%d.fits" % (run, camcol)

        hdulist = fits.open(datafile)
        # nfields = hdulist[0].header['NFIELDS']
        hdu1 = hdulist[1].data

        fieldlist[camcol-1, 0]= np.where(np.max(hdu1[:]['psf_width'])==hdu1[:]['psf_width'])[0][0]
        fieldlist[camcol-1, 1]= np.where(np.min(hdu1[:]['psf_width'])==hdu1[:]['psf_width'])[0][0]
        print('field with max seeing: %d, min seeing: %d'%(fieldlist[camcol-1, 0], fieldlist[camcol-1, 1]))
        for i in range(2):            
            for iBand in range(0, nBand):
                psf = sdsspsf(hdu1, fieldlist[camcol-1, i], iBand, run, camcol)
                psf.fit2vonK_curve_fit(vonK1arcsec)
                if psf.scaleR < -1:
                    psf.fit2vonK_fmin(vonK1arcsec)

                if args.yscale == 'log' or args.yscale == 'logzoom':
                    ax1[iBand, camcol -
                        1].plot(radius * psf.scaleR, np.log10(vonK * psf.scaleV),
                                'b')
                    ax1[iBand, camcol - 1].errorbar(psf.OKprofRadii, psf.OKprofile,
                                                    psf.OKprofileErr, fmt='ok')
                    if args.yscale == 'log':
                        ax1[iBand, camcol - 1].set_xlim(0, 30.0)
                        ax1[iBand, camcol - 1].set_ylim(-6, 0.5)
                        if camcol == 1:
                            text = 'band: %s' % band[iBand]
                            ax1[iBand, camcol -
                                1].text(10, -2, text, fontsize=15, ha='left',
                                        va='center')
                    elif args.yscale == 'logzoom':
                        ax1[iBand, camcol - 1].set_xlim(0, 2)
                        ax1[iBand, camcol - 1].set_ylim(-1.5, 0.1)
                elif args.yscale == 'linear':
                    ax1[iBand, camcol -
                        1].plot(radius * psf.scaleR, vonK * psf.scaleV, color[i])
                    ax1[iBand, camcol - 1].errorbar(psf.OKprofRadii,
                                                    psf.OKprofileLinear,
                                                    psf.OKprofileErrLinear,
                                                    fmt='ok')
                    ax1[iBand, camcol - 1].set_xlim(0, 2.)
                    ax1[iBand, camcol - 1].grid()
                    if camcol == 1:
                        text = 'band: %s' % band[iBand]
                        ax1[iBand, camcol -
                            1].text(0.8, 0.8, text, fontsize=15, ha='left',
                                    va='center')

                if iBand == 0:
                    ax1[iBand, camcol - 1].set_title('camcol=%d' % camcol)
        if camcol ==1:
            titleStr = 'run %d, field %d/%d' % (run, fieldlist[camcol-1, 0], fieldlist[camcol-1, 1])
        else:
            titleStr = '%s, %d/%d' % (titleStr, fieldlist[camcol-1, 0], fieldlist[camcol-1, 1])

    plt.suptitle(titleStr)
    # plt.tight_layout()
    # plt.show()
    plt.savefig('output/run%d_changing_shape_%s.png' % (run, args.yscale))
    end = time.time()
    print('time = %8.2fs' % (end - start))
    sys.exit()


if __name__ == "__main__":
    main()
