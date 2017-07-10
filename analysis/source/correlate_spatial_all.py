import argparse
import numpy as np
from matplotlib import pyplot as plt
import time

from sdssinst import sdssinst

"""
take txt results written out by correlate_spatial.py, combine all the runs
"""

def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_spatial_all.py ---------')
    parser.add_argument('-errorbar', help='if we combine SF from each run equally',
                        action='store_true')
    args = parser.parse_args()
    
    start = time.time()
    
    fwhmStr = 'fwhm'
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    sdss = sdssinst()

    nRow = 2
    nCol = 3
    f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                        sharey='row', figsize=(12, 8))  # a plot is for a run
    # objlist = objlist[objlist[:, 0]<110, :]
    for iBand in range(0, sdss.nBand):
        iRow = np.uint8(np.ceil((iBand+1) / nCol)) - 1
        iCol = np.mod(iBand, nCol)
        for line in objlist:
            run = int(line[0])
            sfname = 'SDSSdata/correlate_spatial/run%d_%s_%s_sf.txt'%(
                run, sdss.band[iBand], fwhmStr)
            txtdata = np.loadtxt(sfname)
            # print('%s'%sfname)
            
            if (run==objlist[0, 0] ):
                mySep = txtdata[0, :]
                nbin = txtdata.shape[1]
                myN = np.zeros(nbin)
                mySF = np.zeros(nbin)
            idx = ~np.isnan( txtdata[1, :])
            mySF[idx] += txtdata[1, idx]**2 * txtdata[2, idx]
            myN[idx] += txtdata[2, idx]

        idx = myN>0
        mySF[idx] = np.sqrt(mySF[idx]/myN[idx])
        mySF[~idx] = -1
            
        if args.errorbar:
            nrun = objlist.shape[0]
            myN[:] = 0 #use this to calculate std around the mean
            for line in objlist:
                run = int(line[0])
                sfname = 'SDSSdata/correlate_spatial/run%d_%s_%s_sf.txt'%(
                    run, sdss.band[iBand], fwhmStr)
                txtdata = np.loadtxt(sfname)
                idx = ~np.isnan( txtdata[1, :])
                myN[idx] += (txtdata[1, idx] - mySF[idx])**2
            myN[idx] = np.sqrt(myN[idx]/nrun)
        
        idx = (mySF>0)
        mySF = mySF[idx]
        mySep = mySep[idx]
        if args.errorbar:
            myN = myN[idx]

        if args.errorbar:
            ax1[iRow, iCol].errorbar(mySep, mySF, myN, fmt='-ok')
        else:
            ax1[iRow, iCol].plot(mySep, mySF, '-ro')
        #ax1[iRow, iCol].set_xlim(0, sdss.nCamcol)
        ax1[iRow, iCol].set_title('runALL, %s, %s' %
                    (fwhmStr, sdss.band[iBand]))
        ax1[iRow, iCol].set_xlabel('Spatial separation')
        #ax1[iRow, iCol].set_ylabel('Covariance (arcsec^2)')
        if iCol == 0:
            ax1[iRow, iCol].set_ylabel('PSF size structure function')
        ax1[iRow, iCol].grid()
        #plt.show()

    for iiCol in range(iCol+1, nCol):
        f.delaxes(ax1[iRow, iiCol])
        
    pngname = 'SDSSdata/correlate_spatial/runALL_%s.png' %(fwhmStr)    
    # plt.tight_layout()
    plt.savefig(pngname)
    plt.close()
        
    end = time.time()
    print('time = %8.2fs' % (end - start))


if __name__ == "__main__":
    main()
