import argparse
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

from sdssinst import sdssinst

"""
for each run (or all runs), each column, draw fwhm(not fwhmvk)
   as a function of wavelength. Overlay all fields in that run.
"""


def main():

    parser = argparse.ArgumentParser(
        description='----- fwhm_lambda.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('yscale', choices=(
        'log', 'linear'), help='yscale of the plots')
    parser.add_argument('-ugri', help='run ugri bands only; w/o z',
                        action='store_true')
    parser.add_argument('-startfield', dest='startfield', default=0, type=int,
                        help='field# to start with')
    parser.add_argument('-endfield', dest='endfield', default=99999, type=int,
                        help='field# to end with (note indexing \
                        starts from 0)')
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
                        action='store_true')
    args = parser.parse_args()

    if args.ugri == 1:
        bands = 'ugri'
    else:
        bands = 'ugriz'
    runNo = args.run
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if runNo > 0:
        # remove all lines but one
        objlist = objlist[objlist[:, 0] == runNo, :]

    sdss = sdssinst()
    runcount = 0

    nRow = 2
    nCol = np.uint8(np.ceil(sdss.nCamcol / nRow))
    f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                          sharey='row', figsize=(12, 8))  # a plot is for a run
    if runNo < 0:
        outTXTfile = 'output/fwhm_lambda/power_%s.txt' % bands
        fidw = open(outTXTfile, 'w')
    xlambda = np.arange(300, 1000, 10)
    for line in objlist:

        run = int(line[0])
        runcount += 1
        print('-- running on run# %d (seq.# %d)---------' % (run, runcount))

        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        idx = (txtdata[:, 0] >= args.startfield) & (
            txtdata[:, 0] < args.endfield)
        txtdata = txtdata[idx, :]
        if args.doubleG:
            fwhmStr = 'fwhm2G'
            fwhm = txtdata[:, 4]
        else:
            fwhmStr = 'fwhm'
            fwhm = txtdata[:, 3] / 1.222  # convert FWHMeff into FWHM
        airmass = txtdata[:, 5]
        fwhm = fwhm/airmass**0.6
        startfield = args.startfield
        endfield = np.uint16(np.max(txtdata[:, 0]))

        if runNo < 0:
            # a plot is for a run
            f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                                  sharey='row', figsize=(12, 8))

        for camcol in range(1, sdss.nCamcol + 1):
            iRow = np.uint8(np.ceil(camcol / nCol)) - 1
            iCol = np.mod(camcol - 1, nCol)
            # print(iRow, iCol)

            # scanning order is r-i-u-z-g
            idx = (txtdata[:, 1] == camcol) & \
                (((txtdata[:, 2] == 0) & (txtdata[:, 0] > startfield + 3) & (
                    txtdata[:, 0] < endfield - 4)) |
                    ((txtdata[:, 2] == 1) & (txtdata[:, 0] < endfield - 8)) |
                 ((txtdata[:, 2] == 2) & (txtdata[:, 0] > startfield + 7)) |
                 ((txtdata[:, 2] == 3) & (txtdata[:, 0] > startfield + 5) & (
                     txtdata[:, 0] < endfield - 2)) |
                 ((txtdata[:, 2] == 4) & (txtdata[:, 0] > startfield + 1) & (
                     txtdata[:, 0] < endfield - 6)))
            if args.ugri == 1:
                idx = idx & (txtdata[:, 2] != 4)  # remove z band data
            fwhmfit = fwhm[idx]
            Leffnparray = np.array(sdss.Leff)
            Lefffit = Leffnparray[np.uint8(txtdata[idx, 2])]
            popt, pcov = optimize.curve_fit(
                lambda Leff, norm: fwhm_lambda_dep(Leff, norm, -0.2),
                Lefffit, fwhmfit, p0=[5])
            norm02 = popt[0]
            y02 = fwhm_lambda_dep(xlambda, norm02, -0.2)
            ax1[iRow, iCol].plot(xlambda, y02, '-b', label='power = -0.2')

            popt, pcov = optimize.curve_fit(
                lambda Leff, norm: fwhm_lambda_dep(Leff, norm, -0.3),
                Lefffit, fwhmfit, p0=[5])
            norm03 = popt[0]
            y03 = fwhm_lambda_dep(xlambda, norm03, -0.3)
            ax1[iRow, iCol].plot(xlambda, y03, '-g', label='power = -0.3')

            popt, pcov = optimize.curve_fit(
                fwhm_lambda_dep,
                Lefffit, fwhmfit, p0=[5, -0.25])
            normfit = popt[0]
            powerfit = popt[1]
            yfit = fwhm_lambda_dep(xlambda, normfit, powerfit)
            ax1[iRow, iCol].plot(xlambda, yfit, '-k',
                                 label='power (fit) = %5.2f' % powerfit)
            print('camcol=%d, power (fit) = %5.2f' % (camcol, powerfit))
            if runNo < 0:
                if camcol == 1:
                    fidw.write('%d\t' % run)
                fidw.write('%5.2f\t' % powerfit)
                if camcol == sdss.nCamcol:
                    fidw.write('\n')

            # scanning order is r-i-u-z-g
            for ufield in range(startfield + 4, endfield - 4 + 1):
                gfield = ufield - 4
                rfield = ufield + 4
                ifield = ufield + 2
                zfield = ufield - 2
                idxu = (txtdata[:, 0] == ufield) & (
                    txtdata[:, 1] == camcol) & (txtdata[:, 2] == 0)
                idxg = (txtdata[:, 0] == gfield) & (
                    txtdata[:, 1] == camcol) & (txtdata[:, 2] == 1)
                idxr = (txtdata[:, 0] == rfield) & (
                    txtdata[:, 1] == camcol) & (txtdata[:, 2] == 2)
                idxi = (txtdata[:, 0] == ifield) & (
                    txtdata[:, 1] == camcol) & (txtdata[:, 2] == 3)
                idxz = (txtdata[:, 0] == zfield) & (
                    txtdata[:, 1] == camcol) & (txtdata[:, 2] == 4)
                if args.ugri == 1:
                    fwhmt = np.hstack(
                        (fwhm[idxu], fwhm[idxg], fwhm[idxr], fwhm[idxi]))
                    ax1[iRow, iCol].plot(sdss.Leff[0:4], fwhmt, '-ro')
                else:
                    fwhmt = np.hstack((fwhm[idxu], fwhm[idxg], fwhm[
                                      idxr], fwhm[idxi], fwhm[idxz]))
                    ax1[iRow, iCol].plot(sdss.Leff, fwhmt, '-ro')
            ax1[iRow, iCol].grid()
            ax1[iRow, iCol].set_title('camcol=%d' % camcol)
            ax1[iRow, iCol].plot(xlambda, y02, '-b')
            ax1[iRow, iCol].plot(xlambda, y03, '-g')
            ax1[iRow, iCol].plot(xlambda, yfit, '-k')

            ax1[iRow, iCol].set_xlabel('Effective wavelength (nm)')
            ax1[iRow, iCol].set_ylabel(fwhmStr)
            ax1[iRow, iCol].legend(loc="upper right", fontsize=10)

            if args.yscale == 'log':
                ax1[iRow, iCol].set_yscale('log')
                ax1[iRow, iCol].set_xscale('log')
                ax1[iRow, iCol].set_xlim([300, 1000])
                ax1[iRow, iCol].set_ylim(
                    [np.min(fwhm) - 0.1, np.max(fwhm) + 0.3])

        if runNo < 0:
            # plt.tight_layout()
            plt.savefig('output/fwhm_lambda/run%d_%s_lambda_%s_%s.png' %
                        (run, fwhmStr, bands, args.yscale))
            plt.close()
    if runNo > 0:
        # plt.tight_layout()
        if (args.startfield == 0 and args.endfield == 99999):
            plt.savefig('output/fwhm_lambda/run%d_%s_lambda_%s_%s.png' %
                        (runNo, fwhmStr, bands, args.yscale))
        else:
            plt.savefig(
                'output/fwhm_lambda/run%d_%s_lambda_%s_%s_%d_%d.png' %
                (runNo, fwhmStr, bands, args.yscale, args.startfield, args.endfield))

        plt.close()

    if runNo < 0:
        fidw.close()


def fwhm_lambda_dep(Leff, norm, power):
    return norm * Leff**power

if __name__ == "__main__":
    main()
