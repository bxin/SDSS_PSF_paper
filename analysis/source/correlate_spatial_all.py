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
    parser.add_argument('-oneplot', help='if we combine SF from all bands',
                        action='store_true')
    args = parser.parse_args()
    
    start = time.time()

    # phosim data
    phosimdata = np.loadtxt('data/phosim_psf_size_sf.txt')
    psr = phosimdata[0, :]
    psatm = phosimdata[1, :]
    psatmE = phosimdata[2, :]
    psall = phosimdata[-2, :]
    psallE = phosimdata[-1, :]
    # CFHT data
    cfhtdata = np.loadtxt('data/CFHT_psf_size_sf.txt')
    cfhtr = cfhtdata[0, :]
    cfhtall = np.mean(cfhtdata[1:,:],0)
    cfhtallE = np.std(cfhtdata[1:,:],0)
    ###
    
    fwhmStr = 'fwhm'
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    sdss = sdssinst()
    
    useRow = 1 # row number in the data files (starts from 0)

    if not args.oneplot:
        nRow = 2
        nCol = 3
        #useRow = 3
        f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                                sharey='row', figsize=(12, 8)) 
        # objlist = objlist[objlist[:, 0]<110, :]
        yymax = 0
        xxmax = 0
    else:
        plt.figure(figsize=(6,4.5))
    
    for iBand in range(0, sdss.nBand):
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
                if args.oneplot:
                    mySFallband = np.zeros(nbin)
                    myErrallband = np.zeros(nbin)
            idx = ~np.isnan( txtdata[useRow, :])
            if (not args.errorbar or args.oneplot):
                mySF[idx] += txtdata[useRow, idx]**2 * txtdata[2, idx]
                myN[idx] += txtdata[2, idx]
            else: #if errorbar is desired for each band:
                mySF[idx] += txtdata[useRow, idx]**2 
                myN[idx] += 1

        idx = myN>0
        mySF[idx] = np.sqrt(mySF[idx]/myN[idx])
        mySF[~idx] = -1
            
        if (args.errorbar or args.oneplot):
            nrun = objlist.shape[0]
            myN[:] = 0 #use this to calculate std around the mean
            for line in objlist:
                run = int(line[0])
                sfname = 'SDSSdata/correlate_spatial/run%d_%s_%s_sf.txt'%(
                    run, sdss.band[iBand], fwhmStr)
                txtdata = np.loadtxt(sfname)
                idx = ~np.isnan( txtdata[useRow, :])
                myN[idx] += (txtdata[useRow, idx] - mySF[idx])**2
            myN[idx] = np.sqrt(myN[idx]/nrun)
        
        idx = (mySF>0)
        mySF = mySF[idx]
        mySep = mySep[idx]
        if (args.errorbar or args.oneplot):
            myN = myN[idx]

        if not args.oneplot:
            iRow = np.uint8(np.ceil((iBand+1) / nCol)) - 1
            iCol = np.mod(iBand, nCol)
            if args.errorbar:
                ax1[iRow, iCol].errorbar(mySep, mySF, myN, fmt='-ok', label='SDSS (%s-band)'%(sdss.band[iBand]))
            else:
                ax1[iRow, iCol].plot(mySep, mySF, '-ok', label='SDSS (%s-band)'%(sdss.band[iBand]))
            #ax1[iRow, iCol].errorbar(psr, psatm, psatmE, fmt='-ob', label='PhoSim/LSST atm only',markersize=3)
            ax1[iRow, iCol].errorbar(psr, psall, psallE, fmt=':xr', label='PhoSim/LSST',markersize=6)
            ax1[iRow, iCol].errorbar(cfhtr, cfhtall, cfhtallE, fmt='--vg', label='CFHT')
            leg = ax1[iRow, iCol].legend(loc="upper left", fontsize=10)
            leg.get_frame().set_alpha(0.5)
            
            #ax1[iRow, iCol].set_xlim(0, sdss.nCamcol)
            #ax1[iRow, iCol].set_title('runALL, %s, %s'%(fwhmStr, sdss.band[iBand]))
            # ax1[iRow, iCol].set_title('%s-band'%(sdss.band[iBand]))
            if iRow == nRow-1: # or iCol == nCol-1:
                ax1[iRow, iCol].set_xlabel('Spatial separation (deg)')
            #ax1[iRow, iCol].set_ylabel('Covariance (arcsec^2)')
            if iCol == 0:
                ax1[iRow, iCol].set_ylabel('PSF size structure function')
            ax1[iRow, iCol].grid()
            #plt.show()
            yymax = np.max(np.hstack((mySF, yymax)))
            xxmax = np.max(np.hstack((mySep, xxmax)))
        else:
            w = 1/myN**2
            mySFallband += mySF**2 * w
            myErrallband += w

    if not args.oneplot:
        for iiCol in range(iCol+1, nCol):
            f.delaxes(ax1[iRow, iiCol])
    
        for iBand in range(0, sdss.nBand):
            iRow = np.uint8(np.ceil((iBand+1) / nCol)) - 1
            iCol = np.mod(iBand, nCol)
            ax1[iRow, iCol].set_ylim(0, yymax*2)
            ax1[iRow, iCol].set_xlim(0,xxmax*1.1)
            
    else:
        mySFallband = mySFallband/myErrallband
        myErrallband = np.sqrt(1/myErrallband)
        plt.errorbar(mySep, mySF, myN, fmt='-ok', label='SDSS')
        plt.errorbar(psr, psall, psallE, fmt=':xr', label='PhoSim/LSST',markersize=6)
        plt.errorbar(cfhtr, cfhtall, cfhtallE, fmt='--vg', label='CFHT')
        leg = plt.legend(loc="upper left", fontsize=10)
        leg.get_frame().set_alpha(0.5)
        
        plt.xlabel('Spatial separation (deg)')
        plt.ylabel('PSF size structure function')
        plt.grid()
            
    pngname = 'SDSSdata/correlate_spatial/runALL_%s.png' %(fwhmStr)        
    # plt.tight_layout()
    plt.savefig(pngname,dpi=500)
    plt.close()
        
    end = time.time()
    print('time = %8.2fs' % (end - start))


if __name__ == "__main__":
    main()
