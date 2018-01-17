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

    if (args.startfield == 0 and args.endfield == 99999):
        fitpname = 'output/correlate_temporal/run%d_%ssf_fitp.txt'%(runNo, fwhmStr)
        pngname = 'output/correlate_temporal/run%d_%ssf_fitp.png'%(runNo, fwhmStr)
    else:
        fitpname = 'output/correlate_temporal/run%d_%s_fld_%d_%d_fitp.txt' %(
            runNo, fwhmStr, args.startfield, args.endfield)
        pngname = 'output/correlate_temporal/run%d_%s_fld_%d_%d_fitp.png' %(
            runNo, fwhmStr, args.startfield, args.endfield)
        
    if args.writefitp:
        fid = open(fitpname, 'w')           
    sdss = sdssinst()
    runcount = 0
    fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8,4))
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
                    dataname = 'output/correlate_temporal/run%d_%ssf_%s_c%d_fitp.txt'%(run, fwhmStr, sdss.band[band], camcol)
                else:
                    dataname = 'output/correlate_temporal/run%d_%ssf_%s_c%d_fld_%d_%d_fitp.txt' %(
                        runNo, fwhmStr, sdss.band[band], camcol, args.startfield, args.endfield)
                    
                bindata = np.loadtxt(dataname)
                try:
                    idx = (bindata[2,:]>1) *(bindata[0, :]<140)
                    if band == 0 and camcol == 1:
                        mySF = np.zeros((np.sum(idx), 5, 6))
                        mySFstd = np.zeros((np.sum(idx), 5, 6))
                
                    mySep = bindata[0, idx]
                    mySF[:, band, camcol-1] = bindata[1, idx]
                    myN = bindata[2, idx]
                    mySFstd[:, band, camcol-1] = bindata[3, idx]
                    
                except IndexError as e:
                    mySep = np.zeros(1)

        c1 =2 #1 # 2
        c2 = 5 #6 # 5
        b2 = 4 #without z-band
        #b2 = 5 #with z-band
        mySFRun = np.mean(mySF[:,:b2,c1-1:c2],(1,2))
        mySFstdRun = np.sqrt(np.mean(mySFstd[:,:b2,c1-1:c2]**2, (1,2))) #approximate, cannot recover original (x_i -x_0) and replace x_0 with mySFRun
        if mySep.shape[0]>2: # we do not fit to <=3 data points
            # mySFstd = mySFstd/100 #playing around, checking why pcov is large
            # curve_fit gives the error
            popt, pcov = optimize.curve_fit(fdt, mySep, mySFRun, p0=[0.1, np.min((30,np.max(mySep))), 0.8],
                                                    bounds = ([0, 5, 0], [0.3, np.max(mySep), 5]),
                                                sigma = mySFstdRun)
            A = popt[0]
            tau = popt[1]
            gamma = popt[2]
            Aerr = np.sqrt(pcov[0, 0])
            tauerr = np.sqrt(pcov[1, 1])
            gammaerr = np.sqrt(pcov[2, 2])
            print('curve_fit: tau = %6.1f\t %6.1f'%(tau, tauerr))
            print('curve_fit: gamma = %6.3f\t %6.3f'%(gamma, gammaerr))
            # fmin does not give the error
            # popt = optimize.fmin(
            #     lambda Atg: fdtChi2(mySep, Atg,mySF, mySFstd), [0.1, 30, 0.8], disp=0)
            # A = popt[0]
            # tau = popt[1]
            # gamma = popt[2]

            # debugging -----
            # fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8,4))
            # ax0.errorbar(mySep, mySFRun, mySFstdRun, fmt = 'ok')
            # myX = np.linspace(0, mySep[0]+mySep[-1], 100)
            # myY = fdt(myX, A, tau, gamma)
            # ax0.plot(myX, myY, 'r-')
            # plt.show()
            if (run == 4874):
                #ax0.errorbar(mySep, mySF, mySFstd, fmt = 'ok')
                ax0.plot(mySep, mySFRun, 'ok')
                myX = np.linspace(0, mySep[0]+mySep[-1], 100)
                myY = fdt(myX, A, tau, gamma)
                chi2 = fdtChi2(mySep, popt, mySFRun, mySFstdRun)
                ax0.plot(myX, myY, 'r-')

                myYDRW = A*np.sqrt(1-(np.exp(-myX/tau)))
                ax0.plot(myX, myYDRW, 'b--')
                # ax0.set_title('run%d, %s-band, camcol=%s' %
                #        (run, sdss.band[band], camcol))
                ax0.grid()
                ax0.set_xlabel(r'$\Delta t$ (minutes)')
                ax0.set_ylabel(r'$<f(\Delta t)>$')
                # plt.show()

            if tau>np.max(mySep):
                tau = 0
                tauerr = 1e8
                gamma = 0
                gammaerr = 1e8                
        else:
            tau = 0
            tauerr = 1e8
            gamma = 0
            gammaerr = 1e8
        # print('fmin: tau = %6.1f'%tau)
        # print('fmin: gamma = %6.3f'%gamma)
        
        print('tau = %6.1f +/- %6.1f, gamma = %6.2f +/- %6.2f\n'% (
            tau, tauerr, gamma, gammaerr))
        if args.writefitp:
            fid.write('%d \t %d \t %6.1f \t %6.1f \t %6.2f \t %6.2f\n'% (
                run, nfields, tau, tauerr, gamma, gammaerr))
        
    if args.writefitp:
        fid.close()

    a = np.loadtxt('output/correlate_temporal/run%d_fwhmsf_fitp.txt' % runNo)
    if runNo<0:
        run = a[:,0]
        totalT = a[:, 1]*36/60
        tau = a[:, 2]
        tauerr = a[:, 3]
        
        idx = totalT>0
        ax1.scatter(totalT[idx], tau[idx], s=10, c='r', edgecolors='r')
        ax1.grid()
        ax1.set_ylabel(r'$\tau$ (minutes)', {'fontsize': 16})
        ax1.set_xlabel('Duration of run (minutes)', {'fontsize': 16})
        ax1.set_xlim(0, 600)
        ax1.set_ylim(0, 120)
        ax1.plot([0,120],[0,120],'--b')
        ax1.plot([0,240],[0,120],'--b')
        
        plt.tight_layout()
        plt.savefig(pngname, dpi=500)
        plt.close()
        
if __name__ == "__main__":
    main()
