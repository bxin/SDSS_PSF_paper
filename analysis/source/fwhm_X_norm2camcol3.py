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
        description='----- fwhm_X.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('yscale', choices=(
        'log', 'linear'), help='yscale of the plots')
    parser.add_argument('-startfield', dest='startfield', default=0, type=int,
                        help='field# to start with')
    parser.add_argument('-endfield', dest='endfield', default=99999, type=int,
                        help='field# to end with (note indexing \
                        starts from 0)')
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
                        action='store_true')
    args = parser.parse_args()

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
        outTXTfile = 'output/fwhm_X/power.txt'
        fidw = open(outTXTfile, 'w')

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
        startfield = args.startfield
        endfield = np.uint16(np.max(txtdata[:, 0]))
        xx = np.arange(min(airmass)-0.01, max(airmass)+0.01, 0.005)
        
        if runNo < 0:
            # a plot is for a run
            f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                                  sharey='row', figsize=(12, 8))

        idxm = txtdata[:, 2]
        fwhmN = np.zeros(txtdata.shape[0])
        for ii in range(txtdata.shape[0]):
            #find the fwhm in (same band/filter, same field number, and camcol = 3)
            # think about a 5x6 array of detectors moving into light; each step move takes 39x2 seconds.
            rx = (txtdata[:, 1] == 3) & (txtdata[:, 0] == txtdata[ii, 0])& (txtdata[:, 2] == txtdata[ii, 2])
            if ((fwhm[ii]/fwhm[rx]==1) & (txtdata[ii, 1] != 3)):
                print('diff col, but happens to have same fwhm: ii=%d, rx=%d'%(ii, np.argmax(rx)))
            fwhmN[ii] = fwhm[ii]/fwhm[rx]
            if ii ==0:
                xnorm = airmass[rx]
                
        for iBand in range(sdss.nBand):
            iRow = np.uint8(np.ceil((iBand+1) / nCol)) - 1
            iCol = np.mod(iBand, nCol)
            # print(iRow, iCol)

            idx = (idxm == iBand) 
            fwhmfit = fwhmN[idx]
            xfit = airmass[idx]

            # get the power=0.6 curve
            y02 = fwhm_X_dep1(xx, xnorm, 0.6)
            ax1[iRow, iCol].plot(xx, y02, '-b', label='power = 0.6')

            # fit for norm and power 
            popt, pcov = optimize.curve_fit(
                lambda x, power:fwhm_X_dep1(x, xnorm, power), 
                xfit, fwhmfit, p0=[0.6])
            powerfit = popt[0]
            yfit = fwhm_X_dep1(xx, xnorm, powerfit)
            ax1[iRow, iCol].plot(xx, yfit, '-k',
                                 label='power (fit) = %5.2f' % powerfit)
            ax1[iRow, iCol].plot(xfit, fwhmfit, 'ro')
            
            print('iBand=%d, power (fit) = %5.2f' % (iBand, powerfit))
            if runNo < 0:
                if iBand == 1:
                    fidw.write('%d\t' % run)
                fidw.write('%5.2f\t' % powerfit)
                if iBand == sdss.nBand:
                    fidw.write('\n')

            fwhmt = np.zeros(sdss.nBand)
            fwhmerr = np.zeros(sdss.nBand)
            for iBand in range(sdss.nBand):
                idxx = (Lefffit == sdss.Leff[iBand])
                aa = fwhm[idx]
                fwhmt[iBand] = np.mean(aa[idxx])
                fwhmerr[iBand] = np.std(aa[idxx])
            ax1[iRow, iCol].errorbar(sdss.Leff, fwhmt, fwhmerr, fmt='ok')
            ax1[iRow, iCol].grid()
            ax1[iRow, iCol].set_title('iBand=%d' % iBand)
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


def fwhm_X_dep(x, norm, power):
    return norm * x**power

def fwhm_X_dep1(x, norm, power):
    return  (x/norm)**power

if __name__ == "__main__":
    main()
