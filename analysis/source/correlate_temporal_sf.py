import argparse
import numpy as np

from astropy.time import Time
from astropy.time import TimeDelta

from sdssinst import sdssinst
import multiprocessing
import time

"""
---for each run, each filter, make plot of cov vs separation
"""

def fdt(dt, A, tau, gamma):
    # print('A = %e, f0 = %f'% (A, f0))
    return A*np.sqrt(1-(np.exp(-(dt/tau)**gamma )))
    #return A*(1-(np.exp(-dt/tau))*gamma )
    # return A*(1-np.exp((-dt/tau)**gamma) )  # (-dt/tau)**gamma -> nan

def fdtChi2(dt, Atg, y, yerr):
    A = Atg[0]
    tau = Atg[1]
    gamma = Atg[2]
    yp = A*np.sqrt(1-(np.exp(-(dt/tau)**gamma)))
    #yp = A*(1-(np.exp(-dt/tau))*gamma)
    #yp = A*(1-np.exp((-dt/tau)**gamma))
    chi2 = np.sum(((yp - y)/yerr)**2)
    # print('A = %e, f0 = %f, chi2 = %f'% (A, f0, chi2))
    return chi2

def sf1run(argList):
    run        =  argList[0] 
    startfield =  argList[1] 
    endfield99 =  argList[2] 
    doubleG    =  argList[3] 
    sdss       =  argList[4] 

    bandlist = np.arange(0,5)
    camcollist = np.arange(1,6+1)
    
    if doubleG:
        fwhmStr = 'fwhm2G'
    else:
        fwhmStr = 'fwhm'

    txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
    txtdata = np.loadtxt(txtfile)
    
    idx = (txtdata[:, 0] >= startfield) & (
        txtdata[:, 0] <= endfield99)
    txtdata = txtdata[idx, :]
    if doubleG:
        fwhm = txtdata[:, 5]
    else:
        fwhm = txtdata[:, 3]  # FWHMeff 
    airmass = txtdata[:, 6]
    mjd = txtdata[:, 7]
    fwhm = fwhm/airmass**0.6

    endfield = np.uint16(np.max(txtdata[:, 0]))
    nfields = endfield - startfield + 1

    for band in bandlist:
        print('band = %d'%band)
        for camcol in camcollist:
            N = nfields

            idx = (txtdata[:, 1] == camcol) & (txtdata[:, 2] == band)
            
            if (startfield == 0 and endfield99 == 99999):
                dataname = 'output/correlate_temporal/run%d_%ssf_%s_c%d_fitp.txt'%(run, fwhmStr, sdss.band[band], camcol)
            else:
                dataname = 'output/correlate_temporal/run%d_%ssf_%s_c%d_fld_%d_%d_fitp.txt' %(
                    run, fwhmStr, sdss.band[band], camcol, startfield, endfield)
                
            fwhmbf = fwhm[idx]
            mjdbf = mjd[idx]
            nn = int(N * (N - 1) /2)
            sfArray = np.zeros((nn, 2)) #0 for separation; 1 for fwhm
            ii = 0
            for i in range(0, N):
                for j in range(i+1, N):
                    t1 = Time(mjdbf[i], format='mjd')
                    t2 = Time(mjdbf[j], format='mjd')
                    dt = t2 - t1
                    sfArray[ii, 0] = abs(TimeDelta(dt).sec)/60
                    if sfArray[ii, 0]>300:
                        break
                    sfArray[ii, 1] = abs(fwhmbf[i] - fwhmbf[j])/(fwhmbf[i] + fwhmbf[j])
                    ii += 1
                    if np.mod(i, 100) == 0 and j==(i+1):
                        print('field=%d, %d'%(i, ii))
            binsize = 10 #10 minutes
            xmin = 0
            xmax = np.ceil(max(sfArray[:, 0])/binsize)*binsize
            nbin = np.int(xmax/binsize)
            myN = np.zeros(nbin) 
            mySF = np.zeros(nbin)
            mySFstd = np.zeros(nbin)
            edge = np.linspace(xmin, xmax, nbin+1)
            mySep = (edge[:-1] + edge[1:])/2
            for i in range(0, nbin):
                idx = (sfArray[:, 0] > edge[i] ) & (sfArray[:, 0]<edge[i+1])
                myN[i] = sum(idx)
                if myN[i] == 0:
                    mySF[i] = np.nan
                    mySFstd[i] = np.nan
                else:
                    mySF[i] = np.mean(sfArray[idx, 1])
                    mySFstd[i] = np.std(sfArray[idx, 1])
            np.savetxt(dataname, np.vstack((mySep, mySF, myN, mySFstd)))
    
def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_temporal.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
                        action='store_true')
    parser.add_argument('-startfield', dest='startfield', default=0, type=int,
                        help='field# to start with')
    parser.add_argument('-endfield', dest='endfield', default=99999, type=int,
                        help='field# to end with (note indexing \
                        starts from 0; endfield is included)')
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
                             sdss))
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
