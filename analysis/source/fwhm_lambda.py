import argparse
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

from sdssinst import sdssinst

# 	for each run (or all runs), each column, draw fwhm(actually fwhmvk)
#			    	 as a function of wavelength. Overlay all fields in that run.

def main():

    parser = argparse.ArgumentParser(
        description='----- fwhm_lambda.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; use -1 for all runs together; -9 for all runs individually')
    args = parser.parse_args()
    
    runNo = args.run
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if runNo>0:
        objlist = objlist[objlist[:,0]==runNo,:] #remove all lines but one
        
    sdss = sdssinst()
    runcount = 0

    nRow = 2
    nCol = np.uint8(np.ceil(sdss.nCamcol/nRow))
    f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                        sharey='row', figsize=(12, 8))  # a plot is for a run

    xlambda = np.arange(300, 1000, 10)
    yfwhm = 0.976*xlambda
    for line in objlist:

        run = int(line[0])
        runcount += 1
        print('-- running on run# %d (seq.# %d)---------' % (run, runcount))

        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        # txtdata = txtdata[txtdata[:, 0]<10, :]
        fwhmvk = txtdata[:, 3]/1.222 #convert FWHMeff into FWHM
        nfields = np.uint16(np.max(txtdata[:, 0]) + 1)
        
        if runNo < -5:
            f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                                sharey='row', figsize=(12, 8))  # a plot is for a run
                        
        for camcol in range(1, sdss.nCamcol + 1):
            iRow = np.uint8(np.ceil(camcol/nCol))-1
            iCol = np.mod(camcol-1, nCol)
            print(iRow, iCol)
            
            idx = (txtdata[:, 1] == camcol)
            fwhmfit = fwhmvk[idx]
            Leffnparray = np.array(sdss.Leff)
            Lefffit = Leffnparray[np.uint8(txtdata[idx, 2])]
            popt, pcov = optimize.curve_fit(
                lambda Leff, norm: fwhm_lambda_dep(-0.2, Leff, norm),
                Lefffit, fwhmfit, p0=[5])
            norm02 = popt[0]
            y02 = fwhm_lambda_dep(-0.2, xlambda, norm02)
            ax1[iRow, iCol].plot(xlambda, y02, '-b')

            popt, pcov = optimize.curve_fit(
                lambda Leff, norm: fwhm_lambda_dep(-0.3, Leff, norm),
                Lefffit, fwhmfit, p0=[5])
            norm03 = popt[0]
            y03 = fwhm_lambda_dep(-0.3, xlambda, norm03)
            ax1[iRow, iCol].plot(xlambda, y03, '-k')
            
            for ifield in range(0, nfields):
                idx = (txtdata[:, 0] == ifield) & (txtdata[:, 1] == camcol)
                ax1[iRow, iCol].plot(sdss.Leff, fwhmvk[idx], '-ro')
            #ax1[iRow, iCol].set_yscale('log')
            #ax1[iRow, iCol].set_xscale('log')
            #ax1[iRow, iCol].set_xlim([300, 1000])
            #ax1[iRow, iCol].set_ylim([1, 3])
            ax1[iRow, iCol].grid()
            ax1[iRow, iCol].set_title('camcol=%d' % camcol)
            ax1[iRow, iCol].plot(xlambda, y02, '-b')
            ax1[iRow, iCol].plot(xlambda, y03, '-k')

            ax1[iRow, iCol].set_xlabel('Effective wavelength (nm)')
            ax1[iRow, iCol].set_ylabel('FWHM')
            
        if runNo < -5:
            plt.tight_layout()
            plt.savefig('output/fwhmvk_lambda/run%d_fwhmvk_lambda.png' % (run))
            plt.close()
    if runNo > -5:        
        # plt.tight_layout()
        if runNo < 0:
            plt.savefig('output/fwhmvk_lambda/runAll_fwhmvk_lambda.png')
        else:
            plt.savefig('output/fwhmvk_lambda/run%d_fwhmvk_lambda.png' % (runNo))

        plt.close()

def fwhm_lambda_dep(power, Leff, norm):
    return norm*Leff**power

if __name__ == "__main__":
    main()
