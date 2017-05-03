
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
        'log', 'linear', 'logzoom'), help='yscale of the plots')
    parser.add_argument('irun', type=int,
                        help='Run Number')
    parser.add_argument('ifield', type=int,
                        help='Field Number')
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

    vdata = np.loadtxt('data/r_vonK_Kolm.txt',
                       unpack='True')
    radius = vdata[0]
    vonK = vdata[1]
    vonK1arcsec = np.vstack((radius, vonK))
    grid1d = np.linspace(-50, 50, 1001)
    vonK2D = np.loadtxt('data/vonK1.0.txt')
    tailPar = np.loadtxt('data/tailPar.txt')

    run = args.irun
    print('------ running on run# %d ---------' % run)
    f, ax1 = plt.subplots(nBand, nCamcol, sharex='col',
                          sharey='row', figsize=(12, 8))  # a plot is for a run

    outdir = "%s/%d/" % (rootName, run)
    for camcol in range(1, nCamcol + 1): #(2,3): #
        print('running on camcol#%d' % camcol)
        datafile = outdir + "photoField-%06d-%d.fits" % (run, camcol)

        hdulist = fits.open(datafile)
        nfields = hdulist[0].header['NFIELDS']
        if args.ifield > nfields:
            print('given field# = %d is too big\n' % args.ifield)
            print('run# %d has %d fields in total\n' % args.irun)
            sys.exit()
        hdu1 = hdulist[1].data

        for iBand in range(0, nBand): # (3,4): #

            psf = sdsspsf(hdu1, args.ifield, iBand, run, camcol)
            psf.tailP = tailPar[iBand*nCamcol + camcol-1]
            
            psf.fit2vonK_curve_fit(vonK1arcsec, vonK2D, grid1d)
            if psf.scaleR < -1:
                psf.fit2vonK_fmin(vonK1arcsec, vonK2D, grid1d)
            psf.fit2convEta_curve_fit(vonK2D, grid1d)
            # psf.fit2conv_curve_fit(vonK2D, grid1d, sigma=0.2)
            # psf.fit2conv_curve_fit(vonK2D, grid1d, tailB = psf.tailB)
            
            print('chi2=%4.1f/%4.1f, chi2lr = %4.1f/%4.1f, chi2hr=%4.1e/%4.1e' %(
                psf.chi2, psf.G2chi2, psf.chi2lr, psf.G2chi2lr, psf.chi2hr, psf.G2chi2hr))

            if args.yscale == 'log' or args.yscale == 'logzoom':
                ax1[iBand, camcol - 1].plot(psf.vR, np.log10(psf.vv), 'b')
                # draw double Gaussian as Red
                # ax1[iBand, camcol -
                #    1].plot(psf.r, psf.LpsfModel + np.log10(psf.scaleV), 'r')
                ax1[iBand, camcol - 1].plot(psf.vvR, np.log10(psf.vvv), 'r')
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
                ax1[iBand, camcol - 1].plot(psf.vR, psf.vv, 'b')
                # draw double Gaussian as Red
                # ax1[iBand, camcol - 1].plot(psf.r,
                #                             psf.psfModel * psf.scaleV, 'r')
                ax1[iBand, camcol - 1].plot(psf.vvR, psf.vvv, 'r')
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

    plt.suptitle('run %d, field %d, Red: vK only, Blue: vK+instrument ' % (run, args.ifield))
    # plt.tight_layout()
    # plt.show()
    plt.savefig('output/run%d_fld%d_psf_vK_2G_%s.png' %
                (run, args.ifield, args.yscale))
    end = time.time()
    print('time = %8.2fs' % (end - start))
    sys.exit()


if __name__ == "__main__":
    main()
