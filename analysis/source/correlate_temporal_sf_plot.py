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
    # bandlist = np.arange(0,4)
    camcollist = np.arange(1,6+1)
    # camcollist = np.arange(2,5+1)
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

        tauAll = np.zeros((5, 6))
        gammaAll = np.zeros((5, 6))
        tauerrAll = np.ones((5, 6))*1e8
        gammaerrAll = np.ones((5, 6))*1e8
                
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
                            # if sfArray[ii, 0]>140:
                            #    break
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
                    bindata = np.loadtxt(dataname)
                    try:
                        mySep = bindata[0, :]
                        mySF = bindata[1, :]
                        myN = bindata[2, :]
                        mySFstd = bindata[3, :]
                        
                        idx = (myN>1) * (mySep<140)
                        mySep = mySep[idx]
                        mySF = mySF[idx]
                        myN = myN[idx]
                        mySFstd = mySFstd[idx]
                
                    except IndexError as e:
                        mySep = np.zeros(1)

                if mySep.shape[0]>2: # we do not fit to <=3 data points
                    # mySFstd = mySFstd/100 #playing around, checking why pcov is large
                    # curve_fit gives the error
                    popt, pcov = optimize.curve_fit(fdt, mySep, mySF, p0=[0.1, np.min((30,np.max(mySep))), 0.8],
                                        bounds = ([0, 5, 0], [0.3, np.max(mySep), 5]),
                                                        sigma = mySFstd)
                    A = popt[0]
                    tau = popt[1]
                    gamma = popt[2]
                    Aerr = np.sqrt(pcov[0, 0])
                    tauerr = np.sqrt(pcov[1, 1])
                    gammaerr = np.sqrt(pcov[2, 2])
                    if tauerr<1e-5:
                        tauerr = 1e8
                    if gammaerr< 1e-5:
                        gammaerr = 1e8

                    # 2-par fits...
                    # popt, pcov = optimize.curve_fit(lambda dt, A, tau:fdt(dt, A, tau, 1), mySep, mySF, p0=[0.1, np.min((30,np.max(mySep)))],
                    #                     bounds = ([0, 5], [0.3, np.max(mySep)]),
                    #                                     sigma = mySFstd)
                    # A = popt[0]
                    # tau = popt[1]
                    # gamma = 1
                    # Aerr = np.sqrt(pcov[0, 0])
                    # tauerr = np.sqrt(pcov[1, 1])
                    # gammaerr = 0
                    # if tauerr<1e-5:
                    #     tauerr = 1e8
                    # if gammaerr< 1e-5:
                    #     gammaerr = 1e8
                        
                    print('curve_fit: tau = %6.1f\t %6.1f'%(tau, tauerr))
                    print('curve_fit: gamma = %6.3f\t %6.3f'%(gamma, gammaerr))
                    # fmin does not give the error
                    # popt = optimize.fmin(
                    #     lambda Atg: fdtChi2(mySep, Atg,mySF, mySFstd), [0.1, 30, 0.8], disp=0)
                    # A = popt[0]
                    # tau = popt[1]
                    # gamma = popt[2]
                    # for debugging -------
                    # fig, (ax0, ax1) = plt.subplots(ncols=2, figsize=(8,4))
                    # ax0.errorbar(mySep, mySF, mySFstd, fmt = 'ok')
                    # myX = np.linspace(0, mySep[0]+mySep[-1], 100)
                    # myY = fdt(myX, A, tau, gamma)
                    # ax0.plot(myX, myY, 'r-')
                    # plt.show()
                    # A = popt[0]
                else:
                    tau = 0
                    tauerr = 1e8
                    gamma = 0
                    gammaerr = 1e8
                # print('fmin: tau = %6.1f'%tau)
                # print('fmin: gamma = %6.3f'%gamma)
                
                tauAll[band,camcol-1] = tau
                gammaAll[band, camcol-1] = gamma
                tauerrAll[band,camcol-1] = tauerr
                gammaerrAll[band, camcol-1] = gammaerr

                if (run == 4874 and band == 2 and camcol == camcollist[0]):
                    #ax0.errorbar(mySep, mySF, mySFstd, fmt = 'ok')
                    ax0.plot(mySep, mySF, 'ok')
                    myX = np.linspace(0, mySep[0]+mySep[-1], 100)
                    myY = fdt(myX, A, tau, gamma)
                    #chi2 = fdtChi2(mySep, popt, mySF, mySFstd)
                    ax0.plot(myX, myY, 'r-')

                    myYDRW = A*np.sqrt(1-(np.exp(-myX/tau)))
                    #ax0.plot(myX, myYDRW, 'b--')
                    # ax0.set_title('run%d, %s-band, camcol=%s' %
                    #        (run, sdss.band[band], camcol))
                    ax0.grid()
                    ax0.set_xlabel(r'$\Delta t$ (minutes)')
                    ax0.set_ylabel(r'$<f(\Delta t)>$')
                    # plt.show()

        w = 1/tauerrAll**2
        tauRun = np.sum(tauAll*w)/np.sum(w)
        tauErrRun =np.sqrt(1/np.sum(w))
        w = 1/gammaerrAll**2
        gammaRun = np.sum(gammaAll*w)/np.sum(w)
        gammaErrRun = np.sqrt(1/np.sum(w))
        print('tau = %6.1f +/- %6.1f, gamma = %6.2f +/- %6.2f\n'% (
            tauRun, tauErrRun, gammaRun, gammaErrRun))
        if args.writefitp:
            fid.write('%d \t %d \t %6.1f \t %6.1f \t %6.2f \t %6.2f\n'% (
                run, nfields, tauRun, tauErrRun, gammaRun, gammaErrRun))
        
    if args.writefitp:
        fid.close()

    a = np.loadtxt('output/correlate_temporal/run%d_fwhmsf_fitp.txt' % runNo)
    if runNo<0:
        run = a[:,0]
        totalT = a[:, 1]*36/60
        tau = a[:, 2]
        tauerr = a[:, 3]
        
        idx = (totalT>0) * (1.5*tau < totalT)

        #plot tau
        ax1.scatter(totalT[idx], tau[idx], s=10, c='r', edgecolors='r')
        ax1.grid()
        ax1.set_ylabel(r'$\tau$ (minutes)', {'fontsize': 16})
        ax1.set_xlabel('Duration of run (minutes)', {'fontsize': 16})
        ax1.set_xlim(0, 600)
        ax1.set_ylim(0, 120)

        # plot gamma
        # gamma = a[:, 4]
        #ax1.scatter(totalT[idx], gamma[idx], s=10, c='r', edgecolors='r')
        # ax1.scatter(tau[idx], gamma[idx], s=10, c='r', edgecolors='r')
        # ax1.grid()
        # ax1.set_ylabel(r'$\gamma$', {'fontsize': 16})
        # ax1.set_xlabel('Duration of run (minutes)', {'fontsize': 16})
        # ax1.set_xlim(0, 600)
        # #ax1.set_ylim(0, 20)
        
        plt.tight_layout()
        plt.savefig(pngname, dpi=500)
        plt.close()
        
if __name__ == "__main__":
    main()
