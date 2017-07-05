import argparse
import numpy as np
from matplotlib import pyplot as plt
#from scipy import optimize
from astropy.coordinates import SkyCoord
import time
import multiprocessing

from sdssinst import sdssinst

"""
---for each run, each filter, make plot of cov vs separation
"""

def func(x, a):
    return a*x+1

def sf1run(argList):
    run        =  argList[0] 
    startfield =  argList[1] 
    endfield   =  argList[2] 
    doubleG    =  argList[3] 
    writesf    =  argList[4] 
    sdss       =  argList[5] 
    
    if doubleG:
        fwhmStr = 'fwhm2G'
    else:
        fwhmStr = 'fwhm'

    txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
    txtdata = np.loadtxt(txtfile)
    idx = (txtdata[:, 0] >= startfield) & (
        txtdata[:, 0] <= endfield)
    txtdata = txtdata[idx, :]
    if doubleG:
        fwhm = txtdata[:, 4]
    else:
        fwhm = txtdata[:, 3]  # FWHMeff 
    airmass = txtdata[:, 5]
    fwhm = fwhm/airmass**0.6
    startfield = startfield
    endfield = np.uint16(np.max(txtdata[:, 0]))
    nfields = endfield - startfield + 1

    nRow = 2
    nCol = 3
    f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                        sharey='row', figsize=(12, 8))  # a plot is for a run
    
    for iBand in range(0, sdss.nBand):
        print('-- running on run# %d , band = %s---------' % (
            run, sdss.band[iBand]))
        if (startfield == 0 and endfield == 99999):
            sfname = 'output/correlate_spatial/run%d_%s_%s_sf.txt'%(
                run, sdss.band[iBand], fwhmStr)
        else:
            sfname = 'output/correlate_spatial/run%d_%s_%s_fld_%d_%d_sf.txt' %(
                run, sdss.band[iBand], fwhmStr, startfield, endfield)
    
        iRow = np.uint8(np.ceil((iBand+1) / nCol)) - 1
        iCol = np.mod(iBand, nCol)

        nn = int(sdss.nCamcol * (sdss.nCamcol - 1) /2)
        sfArray = np.zeros((nfields * nn, 2)) #0 for separation; 1 for fwhm
        ii = 0
        print('iBand = %d, total fields = %d\n'%(iBand, nfields))
        t1 = time.time()
        for ifield in range(0, nfields):
            if (np.mod(ifield, 100) ==0) and (ifield>0):
                t2 = time.time()
                print('ifield=%d, time = %8.2fs'%(ifield, t2 - t1))
                t1 = t2
            idx = (txtdata[:, 2] == iBand) & (txtdata[:, 0] == ifield)
            ra = txtdata[idx, 8]
            dec = txtdata[idx, 9]
            fwhmbf = fwhm[idx] #this is for all combinations of camcol, same band & field
            for i in range(0, sdss.nCamcol):
                for j in range(i+1, sdss.nCamcol):
                    c1 = SkyCoord(ra[i], dec[i], unit = 'deg')
                    c2 = SkyCoord(ra[j], dec[j], unit = 'deg')
                    sfArray[ii, 0] = c1.separation(c2).deg
                    sfArray[ii, 1] = abs(fwhmbf[i] - fwhmbf[j])
                    ii += 1

        nbin = 10
        myN = np.zeros(nbin) 
        mySF = np.zeros(nbin) 
        edge = np.linspace(0, max(sfArray[:, 0]), nbin+1)
        mySep = (edge[:-1] + edge[1:])/2
        for i in range(0, nbin):
            idx = (sfArray[:, 0] > edge[i] ) & (sfArray[:, 0]<edge[i+1])
            myN[i] = sum(idx)
            if myN[i] == 0:
                mySF[i] = np.nan
            else:
                mySF[i] = np.std(sfArray[idx, 1])
        if writesf:
            np.savetxt(sfname, np.vstack((mySep, mySF, myN)))

        idx = ~np.isnan(mySF)
        mySF = mySF[idx]
        mySep = mySep[idx]
        ax1[iRow, iCol].plot(mySep, mySF, '-ro')
        #ax1[iRow, iCol].set_xlim(0, sdss.nCamcol)
        ax1[iRow, iCol].set_title('run%d, %s, %s' %
                    (run, fwhmStr, sdss.band[iBand]))
        ax1[iRow, iCol].set_xlabel('Spatial separation')
        #ax1[iRow, iCol].set_ylabel('Covariance (arcsec^2)')
        ax1[iRow, iCol].set_ylabel('PSF size structure function')
        ax1[iRow, iCol].grid()
        #plt.show()

    if (startfield == 0 and endfield == 99999):
        pngname = 'output/correlate_spatial/run%d_%s.png' %(run, fwhmStr)
    else:
        pngname = 'output/correlate_spatial/run%d_%s_fld_%d_%d.png' %(
            run, fwhmStr, startfield, endfield)

    # plt.tight_layout()
    plt.savefig(pngname)
    plt.close()
    
def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_spatial.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
                        action='store_true')
    parser.add_argument('-writesf', help='write fit parameter',
                        action='store_true')
    parser.add_argument('-startfield', dest='startfield', default=0, type=int,
                        help='field# to start with')
    parser.add_argument('-endfield', dest='endfield', default=99999, type=int,
                        help='field# to end with (note indexing \
                        starts from 0)')
    parser.add_argument('-p', dest='numproc', default=1, type=int,
                    help='Number of Processors to use')                        
    args = parser.parse_args()
    
    start = time.time()

    runNo = args.run
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if runNo > 0:
        # remove all lines but one
        objlist = objlist[objlist[:, 0] == runNo, :]
            
    sdss = sdssinst()
    
    argList = []
    for line in objlist:

        run = int(line[0])
        argList.append((run, args.startfield, args.endfield, args.doubleG,
                            args.writesf, sdss))
        # test
        # sf1run(argList[0])
        
    pool = multiprocessing.Pool(args.numproc)

    pool.map(sf1run, argList)
    pool.close()
    pool.join()

    end = time.time()
    print('time = %8.2fs' % (end - start))


if __name__ == "__main__":
    main()
