import argparse
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

from sdssinst import sdssinst

"""
---for each run, each filter, make plot of cov vs separation
"""


def main():

    parser = argparse.ArgumentParser(
        description='----- correlate_temporal.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; -9 for all runs individually; \
                        makes no sense to plot all runs together')
    parser.add_argument('iBand', type=int, default=0, help='band, 0,1,2,3,4 for ugriz')                        
    parser.add_argument('-doubleG', help='use psf_width from the double Gau fits',
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
            fwhmStr = 'fwhm2G'
            fwhm = txtdata[:, 4]
        else:
            fwhmStr = 'fwhm'
            fwhm = txtdata[:, 3]  # FWHMeff 
        airmass = txtdata[:, 5]
        fwhm = fwhm/airmass**0.6
        startfield = args.startfield
        endfield = np.uint16(np.max(txtdata[:, 0]))
        nfields = endfield - startfield + 1

        nRow = 2
        nCol = np.uint8(np.ceil(sdss.nCamcol / nRow))
        f, ax1 = plt.subplots(nRow, nCol, sharex='col',
                            sharey='row', figsize=(12, 8))  # a plot is for a run
        
        for camcol in range(1, sdss.nCamcol + 1):
            iRow = np.uint8(np.ceil(camcol / nCol)) - 1
            iCol = np.mod(camcol - 1, nCol)

            idx = (txtdata[:, 1] == camcol) & (txtdata[:, 2] == args.iBand)
            
            result = np.correlate(fwhm[idx], fwhm[idx], mode='full')
            autoCorr = result[result.size/2:]
            autoCorr = autoCorr/max(autoCorr)
            ax1[iRow, iCol].plot(autoCorr, marker='.')

            #SFFT = np.fft.fftshift(np.fft.fft(np.fft.fftshift(fwhm[idx])))
            #SFFT = abs(np.fft.fftshift(np.fft.fft(fwhm[idx])))
            #ax1[iRow, iCol].plot(SFFT, marker='.')

            #N = nfields
            ## sample spacing
            #T = 1.0 / N /2
            #x = np.linspace(0.0, N*T, N)
            #y = fwhm[idx]
            #yf = np.fft.fft(y)
            #xf = np.fft.fftfreq(N, T)
            #xf = np.fft.fftshift(xf)
            #yplot = np.fft.fftshift(yf)
            #ax1[iRow, iCol].semilogy(xf, 1.0/N * np.abs(yplot))

            #ax1[iRow, iCol].set_xlim(0, nfields+1)
            ax1[iRow, iCol].set_title('run%d, %s, %s, camcol=%s' %
                        (runNo, fwhmStr, sdss.band[args.iBand], camcol))
            ax1[iRow, iCol].set_xlabel('time (ifield)')
            #ax1[iRow, iCol].set_ylabel('Covariance (arcsec^2)')
            #ax1[iRow, iCol].set_ylabel('Correlation coefficients')
            ax1[iRow, iCol].set_ylabel('FFT')
            ax1[iRow, iCol].grid()
            #plt.show()
            
    # plt.tight_layout()
    if (args.startfield == 0 and args.endfield == 99999):
        plt.savefig('output/correlate_temporal/run%d_%s_%s.png' %
                    (runNo, fwhmStr, sdss.band[args.iBand]))
    else:
        plt.savefig(
            'output/correlate_temporal/run%d_%s_%s_fld_%d_%d.png' %
            (runNo, fwhmStr, sdss.band[args.iBand],args.startfield, args.endfield))

    plt.close()

if __name__ == "__main__":
    main()
