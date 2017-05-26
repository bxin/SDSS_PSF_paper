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

        idxm = np.zeros(txtdata.shape[0])
        fwhmN = np.zeros(txtdata.shape[0])
        for ii in range(txtdata.shape[0]):
            # scanning order is r-i-u-z-g
            # u, normalize using fwhm[r-band]
            if ((txtdata[ii, 2] == 0) & (txtdata[ii, 0] > startfield + 3) & (
                    txtdata[ii, 0] < endfield+1 - 4)):#Nfield = endfield+1
                rx = (txtdata[:, 2] == 2) & (txtdata[:, 0] == txtdata[ii, 0]+4)& (txtdata[:, 1] == txtdata[ii, 1])
                fwhmN[ii] = fwhm[ii]/fwhm[rx]
                idxm[ii]= txtdata[ii, 1]
            # g
            elif ((txtdata[ii, 2] == 1) & (txtdata[ii, 0] < endfield+1 - 8)):
                rx = (txtdata[:, 2] == 2) & (txtdata[:, 0] == txtdata[ii, 0]+8) & (txtdata[:, 1] == txtdata[ii, 1])
                fwhmN[ii] = fwhm[ii]/fwhm[rx]
                idxm[ii]= txtdata[ii, 1]
            # r
            elif ((txtdata[ii, 2] == 2) & (txtdata[ii, 0] > startfield + 7)):
                fwhmN[ii] = 1
                idxm[ii]= txtdata[ii, 1]
            # i
            elif ((txtdata[ii, 2] == 3) & (txtdata[ii, 0] > startfield + 5) & (
                     txtdata[ii, 0] < endfield+1 - 2)):
                rx = (txtdata[:, 2] == 2) & (txtdata[:, 0] == txtdata[ii, 0]+2) & (txtdata[:, 1] == txtdata[ii, 1])
                fwhmN[ii] = fwhm[ii]/fwhm[rx]
                idxm[ii]= txtdata[ii, 1]
            # z
            elif ((txtdata[ii, 2] == 4) & (txtdata[ii, 0] > startfield + 1) & (
                      txtdata[ii, 0] < endfield+1 - 6)):
                rx = (txtdata[:, 2] == 2) & (txtdata[:, 0] == txtdata[ii, 0]+6) & (txtdata[:, 1] == txtdata[ii, 1])
                fwhmN[ii] = fwhm[ii]/fwhm[rx]
                idxm[ii]= txtdata[ii, 1]
                
        for camcol in range(1, sdss.nCamcol + 1):
            iRow = np.uint8(np.ceil(camcol / nCol)) - 1
            iCol = np.mod(camcol - 1, nCol)
            # print(iRow, iCol)

            idx = (idxm == camcol) 
            if args.ugri == 1:
                idx = idx & (txtdata[:, 2] != 4)  # remove z band data
                sdss.nBand = 4
                sdss.Leff = sdss.Leff[:4]
            fwhmfit = fwhmN[idx]
            Leffnparray = np.array(sdss.Leff)
            Lefffit = Leffnparray[np.uint8(txtdata[idx, 2])]
            y02 = fwhm_lambda_dep1(xlambda, -0.2)
            ax1[iRow, iCol].plot(xlambda, y02, '-b', label='power = -0.2')

            y03 = fwhm_lambda_dep1(xlambda, -0.3)
            ax1[iRow, iCol].plot(xlambda, y03, '-g', label='power = -0.3')

            popt, pcov = optimize.curve_fit(
                fwhm_lambda_dep1,
                Lefffit, fwhmfit, p0=[-0.25])
            powerfit = popt[0]
            yfit = fwhm_lambda_dep1(xlambda, powerfit)
            ax1[iRow, iCol].plot(xlambda, yfit, '-k',
                                 label='power (fit) = %5.2f' % powerfit)
            print('camcol=%d, power (fit) = %5.2f' % (camcol, powerfit))
            if runNo < 0:
                if camcol == 1:
                    fidw.write('%d\t' % run)
                fidw.write('%5.2f\t' % powerfit)
                if camcol == sdss.nCamcol:
                    fidw.write('\n')

            fwhmt = np.zeros(sdss.nBand)
            fwhmerr = np.zeros(sdss.nBand)
            for iBand in range(sdss.nBand):
                idxx = (Lefffit == sdss.Leff[iBand])
                aa = fwhmN[idx]
                fwhmt[iBand] = np.mean(aa[idxx])
                fwhmerr[iBand] = np.std(aa[idxx])
            ax1[iRow, iCol].errorbar(sdss.Leff, fwhmt, fwhmerr, fmt='ok')
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
                    [np.min(fwhmN) - 0.1, np.max(fwhmN) + 0.3])

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

def fwhm_lambda_dep1(Leff, power):
    return  (Leff/616.6)**power

if __name__ == "__main__":
    main()
