import argparse
import numpy as np
from matplotlib import pyplot as plt
# from matplotlib import mlab
from scipy import optimize
from scipy import fftpack
from astropy.time import Time
from astropy.time import TimeDelta

from sdssinst import sdssinst
from astroML.fourier import PSD_continuous

"""
---for each run, each filter, make plot of cov vs separation
"""

def fdt(dt, A, tau, gamma):
    # print('A = %e, f0 = %f'% (A, f0))
    return A*(1-(np.exp(-dt/tau))**gamma )
    #return A*(1-(np.exp(-dt/tau))*gamma )
    # return A*(1-np.exp((-dt/tau)**gamma) )  # (-dt/tau)**gamma -> nan

def fdtChi2(dt, Atg, y, yerr):
    A = Atg[0]
    tau = Atg[1]
    gamma = Atg[2]
    yp = A*(1-(np.exp(-dt/tau))**gamma)
    #yp = A*(1-(np.exp(-dt/tau))*gamma)
    #yp = A*(1-np.exp((-dt/tau)**gamma))
    chi2 = np.sum(((yp - y)/yerr)**2)
    # print('A = %e, f0 = %f, chi2 = %f'% (A, f0, chi2))
    return chi2

def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_temporal.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('iBand', type=int, default=0, help='band, 0,1,2,3,4 for ugriz,\
                        -9 for looping over all bands')
    parser.add_argument('iCamcol', type=int, default=0, help='camcol, 1-6,\
                        -9 for looping over all camcol')
    parser.add_argument('-writesf', help='write sf data',
                        action='store_true')
    parser.add_argument('-writefitp', help='write fit parameter',
                        action='store_true')
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
                        action='store_true')
    parser.add_argument('-startfield', dest='startfield', default=0, type=int,
                        help='field# to start with')
    parser.add_argument('-endfield', dest='endfield', default=99999, type=int,
                        help='field# to end with (note indexing \
                        starts from 0; endfield is included)')
    args = parser.parse_args()

    runNo = args.run
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    bandlist = np.arange(0,5)
    camcollist = np.arange(1,6+1)
    if runNo > 0:
        # remove all lines but one
        objlist = objlist[objlist[:, 0] == runNo, :]
    if args.iBand > -1:
        bandlist = bandlist[bandlist == args.iBand]
    if args.iCamcol>-1:
        camcollist = camcollist[camcollist == args.iCamcol]
        
    if args.doubleG:
        fwhmStr = 'fwhm2G'
    else:
        fwhmStr = 'fwhm'

    sdss = sdssinst()
    runcount = 0
    for line in objlist:

        run = int(line[0])
        runcount += 1
        print('-- running on run# %d (seq.# %d)---------' % (run, runcount))

        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        if np.mod(args.endfield - args.startfield + 1, 2)==0:
            idx = (txtdata[:, 0] >= args.startfield) & (
                txtdata[:, 0] <= args.endfield)
        else:
            idx = (txtdata[:, 0] >= args.startfield) & (
                txtdata[:, 0] <= args.endfield-1)
        txtdata = txtdata[idx, :]
        if args.doubleG:
            fwhm = txtdata[:, 5]
        else:
            fwhm = txtdata[:, 3]  # FWHMeff 
        airmass = txtdata[:, 6]
        mjd = txtdata[:, 7]
        fwhm = fwhm/airmass**0.6
        startfield = args.startfield
        endfield = np.uint16(np.max(txtdata[:, 0]))
        nfields = endfield - startfield + 1

        for band in bandlist:
            print('band = %d'%band)
            for camcol in camcollist:
                N = nfields
    
                idx = (txtdata[:, 1] == camcol) & (txtdata[:, 2] == band)
                
                if (args.startfield == 0 and args.endfield == 99999):
                    dataname = 'output/correlate_temporal/run%d_%ssf_%s_c%d_fitp.txt'%(runNo, fwhmStr, sdss.band[band], camcol)
                    fitpname = 'output/correlate_temporal/run%d_%s_fitp.txt'%(runNo, fwhmStr)
                else:
                    dataname = 'output/correlate_temporal/run%d_%ssf_%s_c%d_fld_%d_%d_fitp.txt' %(
                        runNo, fwhmStr, sdss.band[band], camcol, args.startfield, args.endfield)
                    fitpname = 'output/correlate_temporal/run%d_%s_fld_%d_%d_fitp.txt' %(
                        runNo, fwhmStr, args.startfield, args.endfield)
                    
                if args.writesf:
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
                            if sfArray[ii, 0]>140:
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
                else:
                    txtdata = np.loadtxt(dataname)
                    mySep = txtdata[0, :]
                    mySF = txtdata[1, :]
                    myN = txtdata[2, :]
                    mySFstd = txtdata[3, :]

                # mySFstd = mySFstd/100 #playing around, checking why pcov is large
                # curve_fit gives the error
                popt, pcov = optimize.curve_fit(fdt, mySep, mySF, p0=[0.1, 30, 0.8],
                                                        bounds = ([0, 5, 0], [0.3, 300, 5]),
                                                    sigma = mySFstd)
                A = popt[0]
                tau = popt[1]
                gamma = popt[2]
                Aerr = np.sqrt(pcov[0, 0])
                tauerr = np.sqrt(pcov[1, 1])
                gammaerr = np.sqrt(pcov[2, 2])
                print('curve_fit: tau = %6.1f\t %6.1f'%(tau, tauerr))
                print('curve_fit: gamma = %6.3f\t %6.3f'%(gamma, gammaerr))
                # fmin does not give the error
                popt = optimize.fmin(
                    lambda Atg: fdtChi2(mySep, Atg,mySF, mySFstd), [0.1, 30, 0.8], disp=0)
                A = popt[0]
                tau = popt[1]
                gamma = popt[2]
                print('fmin: tau = %6.1f'%tau)
                print('fmin: gamma = %6.3f'%gamma)
                
                if args.writefitp:
                    fid = open(fitpname, 'w')     
                    fid.write('%5.3f\t %5.3f\n' % (A, Aerr))
                    fid.write('%5.0f\t %5.0f\n' % (tau, tauerr))
                    fid.write('%5.2f\t %5.2f\n' % (gamma, gammaerr))
                    
                plt.figure(figsize=(6,4.5))
                plt.errorbar(mySep, mySF, mySFstd, fmt = 'ok')
                myX = np.linspace(0, mySep[0]+mySep[-1], 100)
                myY = fdt(myX, A, tau, gamma)
                chi2 = fdtChi2(mySep, popt, mySF, mySFstd)
                plt.plot(myX, myY, 'r-')
                
                plt.title('run%d, %s, %s, camcol=%s' %
                            (run, fwhmStr, sdss.band[band], camcol))
                plt.grid()
                # plt.show()
                
            if (args.startfield == 0 and args.endfield == 99999):
                pngname = 'output/correlate_temporal/run%d_%s_sf_%s.png' %(
                    run, fwhmStr, sdss.band[band])
            else:
                pngname = 'output/correlate_temporal/run%d_%s_sf_%s_fld_%d_%d.png' %(
                    run, fwhmStr, sdss.band[band], args.startfield, args.endfield)
                
            # plt.tight_layout()
            plt.savefig(pngname)
            plt.close()
        
    if args.writefitp:
        fid.close()
        
if __name__ == "__main__":
    main()
