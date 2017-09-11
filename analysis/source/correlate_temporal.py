import argparse
import numpy as np
from matplotlib import pyplot as plt
# from matplotlib import mlab
from scipy import optimize
from scipy import fftpack

from sdssinst import sdssinst
from astroML.fourier import PSD_continuous

"""
---for each run, each filter, make plot of cov vs separation
"""

def logpsdfunc(f, A, f0):
    # print('A = %e, f0 = %f'% (A, f0))
    return np.log10(A/(1+(f/f0)**2))

def logpsdfuncChi2(f, Af0, y):
    A = Af0[0]
    f0 = Af0[1]
    yp = np.log10(A/(1+(f/f0)**2))
    chi2 = np.sum((yp - y)**2)
    # print('A = %e, f0 = %f, chi2 = %f'% (A, f0, chi2))
    return chi2

def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_temporal.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('iBand', type=int, default=0, help='band, 0,1,2,3,4 for ugriz\
                        -9 for looping over all bands')
    parser.add_argument('-type', dest='type',
                        choices=('psd', 'autocor'), default = 'psd', 
                        help='calculate psd or autocorrelation')
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
    if runNo > 0:
        # remove all lines but one
        objlist = objlist[objlist[:, 0] == runNo, :]
    if args.iBand > -1:
        bandlist = bandlist[bandlist == args.iBand]
        
    if args.doubleG:
        fwhmStr = 'fwhm2G'
    else:
        fwhmStr = 'fwhm'

    if (args.startfield == 0 and args.endfield == 99999):
        fitpname = 'output/correlate_temporal/run%d_%s_fitp.txt'%(runNo, fwhmStr)
    else:
        fitpname = 'output/correlate_temporal/run%d_%s_fld_%d_%d_fitp.txt' %(
            runNo, fwhmStr, args.startfield, args.endfield)

    txtdata= np.loadtxt('data/opsim_seeing_psd.txt')
    opsimf = txtdata[:, 0]
    opsimpsd = txtdata[:, 1]
    
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
        if np.mod(args.endfield - args.startfield + 1, 2)==0:
            idx = (txtdata[:, 0] >= args.startfield) & (
                txtdata[:, 0] <= args.endfield)
        else:
            idx = (txtdata[:, 0] >= args.startfield) & (
                txtdata[:, 0] <= args.endfield-1)
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
        nCol = np.uint8(np.ceil(sdss.nCamcol / nRow))
        for band in bandlist:
            print('band = %d'%band)
            f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                            sharey='row', figsize=(12, 8))  # a plot is for a run
            for camcol in range(1, sdss.nCamcol + 1):
                iRow = np.uint8(np.ceil(camcol / nCol)) - 1
                iCol = np.mod(camcol - 1, nCol)
                N = nfields
                dt = 36
    
                idx = (txtdata[:, 1] == camcol) & (txtdata[:, 2] == band)
    
                if args.type == 'psd':
                    #SFFT = np.fft.fftshift(np.fft.fft(np.fft.fftshift(fwhm[idx])))
                    # the line above gives same results as below, except for the abs(), i.e., the
                    # fftshift immediately outside of fwhm doesn't do anything
                    # SFFT = abs(np.fft.fftshift(np.fft.fft(fwhm[idx])))
                    # Neither is needed if we do the below.
                    
                    ## ax1[iRow, iCol].plot(SFFT, marker='.')
                    
                    # 4/23/2017, this part produces same results as the section below (when dt=1)
                    # T = 1.0 / N /2
                    # x = np.linspace(0.0, N*T, N)
                    # y = fwhm[idx]
                    # yf = np.fft.fft(y)
                    # xf = np.fft.fftfreq(N, T)
                    # xf = np.fft.fftshift(xf)
                    # yplot = np.fft.fftshift(yf)
                    # ax1[iRow, iCol].semilogy(xf, 1.0/N * np.abs(yplot))
                    
                    # based on code from
                    # http://www.astroml.org/book_figures/chapter10/fig_LIGO_power_spectrum.html
                    # compute PSD using simple FFT
                    df = 1./ (N * dt)
                    aa = abs(fftpack.fft(fwhm[idx]))**2 * dt/N
                    PSD = aa[:int(N / 2)]
                    PSD += aa[-1:-int(N / 2)-1:-1]
                    f = df * np.arange(int(N / 2))

                    # npsd = 64 # 2048 # 16
                    # # compute PSD using Welch's method -- no window function
                    # PSDW1, fW1 = mlab.psd(fwhm[idx], NFFT=npsd, Fs=1. / dt,
                    #           window=mlab.window_none, noverlap=npsd/2)
                    # dfW1 = fW1[1] - fW1[0]
                    # 
                    # # compute PSD using Welch's method -- hanning window function
                    # PSDW2, fW2 = mlab.psd(fwhm[idx], NFFT=npsd, Fs=1. / dt,
                    #                           window=mlab.window_hanning, noverlap=npsd/2)
                    # dfW2 = fW2[1] - fW2[0]
                    cutoff = f>0
                    f = f[cutoff]
                    PSD = PSD[cutoff]

                    popt, pcov = optimize.curve_fit(logpsdfunc, f, np.log10(PSD), p0=[1e5, 1e-3],
                                                        bounds = ([1, 1e-4], [1e7, 1e-2]))
                    # popt = optimize.fmin(lambda Af0: logpsdfuncChi2(f, Af0, np.log10(PSD)), [1e5, 1e-3], disp=1)
                                                   
                    # sigma=myErr, absolute_sigma=True)
                    tau = 1/(2*np.pi*popt[1])
                    # print('A=%e, f0=%e, tau = %5.0f s' %(popt[0], popt[1], tau))
                    if args.writefitp:
                        fid.write('%d\t %d\t %d\t %e\t %e\t %6.1f \n'%(run, band, camcol, popt[0], popt[1], tau/60))

                    myX = np.hstack((f[0]/2, f, f[-1]*2))
                    myY = 10**logpsdfunc(myX, popt[0], popt[1])
                    chi2 = logpsdfuncChi2(f, popt, np.log10(PSD))
                    # popt = [1e4, 6e-4]
                    # myY1 = 10**logpsdfunc(myX, popt[0], popt[1])
                    ax1[iRow, iCol].loglog(opsimf*3600, opsimpsd, linestyle = 'None', marker='.', color='y', markersize=10)#, c='#AAAAAA')
                    ax1[iRow, iCol].loglog(f*3600, PSD, linestyle = 'None', marker='.', color='k', markersize=10)#, c='#AAAAAA')
                    # ax1[iRow, iCol].loglog(fW1, PSDW1,'-k')
                    # ax1[iRow, iCol].loglog(fW2, PSDW2,'-r')
                    ax1[iRow, iCol].loglog(myX*3600, myY, 'r-')
                    # ax1[iRow, iCol].loglog(myX, myY1, 'b-')
                    
                    ax1[iRow, iCol].set_xlim(min(myX*3600), max(myX*3600))
                    # ax1[iRow, iCol].set_xlabel('Frequency (Hz)')
                    ax1[iRow, iCol].set_xlabel('Frequency (1/hour)')
                    ax1[iRow, iCol].set_ylabel('PSD (arcsec$^2$ second)')

                    print('%d\t %d\t %d\t %e\t %e\t %6.1f \n'%(run, band, camcol, popt[0], popt[1], tau/60))
                    
                    # below, use PSD_continuous and see how different things are
                    # hj_N = fwhm[idx]
                    # if np.mod(N, 2)==1:
                    #     N -= 1
                    #     hj_N = hj_N[:-1]
                    # tj = np.linspace(-N/2*dt, N/2*dt, N)
                    # fc, PSDc = PSD_continuous(tj, hj_N)
                    # cutoff = fc>0
                    # fc = fc[cutoff]
                    # PSDc = PSDc[cutoff]/dt/N
                    # popt, pcov = optimize.curve_fit(logpsdfunc, fc, np.log10(PSDc), p0=[1e5, 1e-3],
                    #                                     bounds = ([1, 1e-5], [1e8, 1e-2]))
                    # # popt = optimize.fmin(lambda Af0: logpsdfuncChi2(f, Af0, np.log10(PSD)), [1e5, 1e-3], disp=1)
                    #                                
                    # # sigma=myErr, absolute_sigma=True)
                    # tau = 1/(2*np.pi*popt[1])
                    # # print('A=%e, f0=%e, tau = %5.0f s' %(popt[0], popt[1], tau))
                    # if args.writefitp:
                    #     fid.write('c: %d\t %d\t %d\t %e\t %e\t %6.1f \n'%(run, band, camcol, popt[0], popt[1], tau/60))
                    # print('c: %d\t %d\t %d\t %e\t %e\t %6.1f \n'%(run, band, camcol, popt[0], popt[1], tau/60))
                    # myX = np.hstack((f[0]/2, f, f[-1]*2))
                    # myY = 10**logpsdfunc(myX, popt[0], popt[1])
                    # ax1[iRow, iCol].loglog(fc, PSDc, linestyle = 'None', marker='.', color='g', markersize=10)#, c='#AAAAAA')
                    # ax1[iRow, iCol].loglog(myX*3600, myY, 'b-')
                    # 
                    # chi21 = logpsdfuncChi2(f, popt, np.log10(PSD))

                elif args.type == 'autocor':
                    result = np.correlate(fwhm[idx], fwhm[idx], mode='full')
                    autoCorr = result[result.size/2:]
                    autoCorr = autoCorr/max(autoCorr)
                    x = np.arange(N) * dt
                    ax1[iRow, iCol].plot(x, autoCorr, marker='.')
    
                    ax1[iRow, iCol].set_xlabel('Time (seconds)')
                    ax1[iRow, iCol].set_ylabel('Correlation coefficients')
                
                ax1[iRow, iCol].set_title('run%d, %s, %s, camcol=%s' %
                            (run, fwhmStr, sdss.band[band], camcol))
                ax1[iRow, iCol].grid()
                #plt.show()
                
            if (args.startfield == 0 and args.endfield == 99999):
                pngname = 'output/correlate_temporal/run%d_%s_%s_%s.png' %(
                    run, fwhmStr, args.type, sdss.band[band])
            else:
                pngname = 'output/correlate_temporal/run%d_%s_%s_%s_fld_%d_%d.png' %(
                    run, fwhmStr, args.type, sdss.band[band], args.startfield, args.endfield)
                
            # plt.tight_layout()
            plt.savefig(pngname)
            plt.close()
        
    if args.writefitp:
        fid.close()
        
if __name__ == "__main__":
    main()
