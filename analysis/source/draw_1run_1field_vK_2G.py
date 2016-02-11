
import sys
import time
import argparse

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

from sdsspsf import sdsspsf


def main():
    parser = argparse.ArgumentParser(
        description='----- draw_1run_1field_vK_2G.py ---------')
    parser.add_argument('yscale', choices=(
        'log', 'linear'), help='yscale of the plots')
    parser.add_argument('-run', dest='run', type=int, default=94,
                        help='run number we want to look at default=94')
    parser.add_argument('-field', dest='field', type=int, default=0,
                        help='field number we want to look at default=0')
    args = parser.parse_args()

    start = time.time()
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if args.run not in objlist[0]:
        print('run# %d is not in Stripe82RunList.dat\n' % args.run)
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

    vdata = np.loadtxt('data/r_vonK_Kolm.txt',
                       unpack='True')
    radius = vdata[0]
    vonK = vdata[1]
    vonK1arcsec = np.vstack((radius, vonK))

    run = args.run
    print('------ running on run# %d ---------' % run)
    f, ax1 = plt.subplots(nBand, nCamcol, sharex='col',
                          sharey='row', figsize=(12, 8))  # a plot is for a run

    outdir = "%s/%d/" % (rootName, run)
    for camcol in range(1, nCamcol + 1):
        print('running on camcol#%d' % camcol)
        datafile = outdir + "photoField-%06d-%d.fits" % (run, camcol)

        hdulist = fits.open(datafile)
        nfields = hdulist[0].header['NFIELDS']
        if args.field > nfields:
            print('given field# = %d is too big\n' % args.field)
            print('run# %d has %d fields in total\n' % args.run)
            sys.exit()
        hdu1 = hdulist[1].data

        for iBand in range(0, nBand):

            psf = sdsspsf(hdu1, args.field, iBand, run, camcol)
            psf.fit2vonK_curve_fit(vonK1arcsec)
            # psf.fit2vonK_fminbound(vonK1arcsec)

            if args.yscale == 'log':
                ax1[iBand, camcol -
                    1].plot(radius * psf.scaleR, np.log10(vonK * psf.scaleV),
                            'b')
                ax1[iBand, camcol -
                    1].plot(psf.r, psf.LpsfModel + np.log10(psf.scaleV), 'r')
                ax1[iBand, camcol - 1].errorbar(psf.OKprofRadii, psf.OKprofile,
                                                psf.OKprofileErr, fmt='ok')
                ax1[iBand, camcol - 1].set_xlim(0, 30.0)
                ax1[iBand, camcol - 1].set_ylim(-6, 0.5)
                if camcol == 1:
                    text = 'band: %s' % band[iBand]
                    ax1[iBand, camcol -
                        1].text(10, -2, text, fontsize=15, ha='left',
                                va='center')
            elif args.yscale == 'linear':
                ax1[iBand, camcol -
                    1].plot(radius * psf.scaleR, vonK * psf.scaleV, 'b')
                ax1[iBand, camcol - 1].plot(psf.r,
                                            psf.psfModel * psf.scaleV, 'r')
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

    plt.suptitle('run %d, field %d, blue: vK, Red: 2G ' % (run, args.field))
    # plt.tight_layout()
    # plt.show()
    plt.savefig('output/run%dpsf_vK_2G_%s.png' % (run, args.yscale))
    end = time.time()
    print('time = %8.2fs' % (end - start))
    sys.exit()


if __name__ == "__main__":
    main()
