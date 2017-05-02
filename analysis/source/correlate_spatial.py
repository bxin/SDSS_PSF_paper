import argparse
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

from sdssinst import sdssinst

"""
---for each run, each filter, make plot of cov vs separation
"""

def func(x, a):
    return a*x+1

def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_spatial.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
                        action='store_true')
    parser.add_argument('-writefitp', help='write fit parameter',
                        action='store_true')
    parser.add_argument('-startfield', dest='startfield', default=0, type=int,
                        help='field# to start with')
    parser.add_argument('-endfield', dest='endfield', default=99999, type=int,
                        help='field# to end with (note indexing \
                        starts from 0)')
    args = parser.parse_args()

    runNo = args.run
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if runNo > 0:
        # remove all lines but one
        objlist = objlist[objlist[:, 0] == runNo, :]
    if args.doubleG:
        fwhmStr = 'fwhm2G'
    else:
        fwhmStr = 'fwhm'
            
    if (args.startfield == 0 and args.endfield == 99999):
        fitpname = 'output/correlate_spatial/run%d_%s_fitp.txt'%(runNo, fwhmStr)
    else:
        fitpname = 'output/correlate_spatial/run%d_%s_fld_%d_%d_fitp.txt' %(
            runNo, fwhmStr, args.startfield, args.endfield)
        
    if args.writefitp:
        fid = open(fitpname, 'w')        
    sdss = sdssinst()
    runcount = 0
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
            fwhm = txtdata[:, 4]
        else:
            fwhm = txtdata[:, 3]  # FWHMeff 
        airmass = txtdata[:, 5]
        fwhm = fwhm/airmass**0.6
        startfield = args.startfield
        endfield = np.uint16(np.max(txtdata[:, 0]))
        nfields = endfield - startfield + 1

        nRow = 2
        nCol = 3
        f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                            sharey='row', figsize=(12, 8))  # a plot is for a run
        
        for iBand in range(0, sdss.nBand):
            iRow = np.uint8(np.ceil((iBand+1) / nCol)) - 1
            iCol = np.mod(iBand, nCol)

            fwhmArray = np.zeros((nfields, sdss.nCamcol))
            for camcol in range(1, sdss.nCamcol + 1):
                idx = (txtdata[:, 1] == camcol) & (txtdata[:, 2] == iBand)
                fwhmArray[:, camcol-1] = fwhm[idx]
            #myCovSquare=np.cov(fwhmArray, rowvar=0)
            myCorrCoef=np.corrcoef(fwhmArray, rowvar=0)
            myCor = np.zeros(15) #15 = 1+2+3+4+5
            mySep = np.zeros(15)
            myErr = np.zeros(15)
            i = 0
            for iCamcol in range(1, sdss.nCamcol + 1):
                for jCamcol in range(iCamcol+1, sdss.nCamcol + 1):
                    myCor[i] = myCorrCoef[iCamcol-1, jCamcol-1]
                    mySep[i] = jCamcol - iCamcol
                    # https://www.mathworks.com/matlabcentral/newsreader/view_thread/107045
                    myErr[i] = (1-myCor[i]**2)/np.sqrt(nfields-1)
                    i += 1

            popt, pcov = optimize.curve_fit(func, mySep, myCor, p0=[-0.005],
                sigma=myErr, absolute_sigma=True)
            myX = np.linspace(0, sdss.nCamcol, 100)
            myY = func(myX, popt[0])
            if args.writefitp:
                fid.write('%d\t %d\t %6.4f \n'%(run, iBand, popt[0]))
                
            wSep = np.zeros(sdss.nCamcol-1)
            wCor = np.zeros(sdss.nCamcol-1)
            wErr = np.zeros(sdss.nCamcol-1)
            for iSep in range(1, sdss.nCamcol):
                myidx = mySep == iSep
                wSep[iSep-1] = iSep
                wCor[iSep-1] = sum(myCor[myidx]/myErr[myidx]**2)/sum(1/myErr[myidx]**2)
                wErr[iSep-1] = np.sqrt(1/sum(1/myErr[myidx]**2))
            
            ax1[iRow, iCol].errorbar(wSep, wCor, wErr, fmt = 'ok') #,markersize=15)
            ax1[iRow, iCol].plot(myX, myY, '-r')
            ax1[iRow, iCol].set_xlim(0, sdss.nCamcol)
            ax1[iRow, iCol].set_title('run%d, %s, %s' %
                        (run, fwhmStr, sdss.band[iBand]))
            ax1[iRow, iCol].set_xlabel('Spatial separation')
            #ax1[iRow, iCol].set_ylabel('Covariance (arcsec^2)')
            ax1[iRow, iCol].set_ylabel('Correlation coefficients')
            ax1[iRow, iCol].grid()
            #plt.show()

        if (args.startfield == 0 and args.endfield == 99999):
            pngname = 'output/correlate_spatial/run%d_%s.png' %(run, fwhmStr)
        else:
            pngname = 'output/correlate_spatial/run%d_%s_fld_%d_%d.png' %(
                run, fwhmStr, args.startfield, args.endfield)

        # plt.tight_layout()
        plt.savefig(pngname)
        plt.close()
        
    if args.writefitp:
        fid.close()
        
if __name__ == "__main__":
    main()
