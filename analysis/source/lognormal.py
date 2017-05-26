import time
import argparse

import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

def lognormal(x, mu, sigma, A):

    return A/sigma/x/np.sqrt(2*np.pi)*np.exp(-(np.log(x)-mu)**2/2/sigma**2)

def main():
    parser = argparse.ArgumentParser(
        description='----- lognormal.py ---------')    
    parser.add_argument('iBand', type=int, default=0, help='band, 0,1,2,3,4 for ugriz\
                        -9 for looping over all bands')
    parser.add_argument('-readhist', help='read histogram data from txt file',
                        action='store_true')
    args = parser.parse_args()
    
    start = time.time()
    bandL = {}
    bandL[0] = "u"
    bandL[1] = "g"
    bandL[2] = "r"
    bandL[3] = "i"
    bandL[4] = "z"

    totalFields = 31580 * 6
    bw = 2.7*0.8/totalFields**(1./3)
    nbins = int(4 / bw)
    runcount = 0
    
    objlist = np.loadtxt('data/Stripe82RunList.dat')

    histtxt = 'output/lognormalHist_%d.txt' % args.iBand
    if args.readhist:
        a = np.loadtxt(histtxt)
        x = a[:, 0]
        y = a[:, 1]
    else:
        seeing = np.zeros(totalFields)
        ii = 0
        for line in objlist:    
            run = int(line[0])
            runcount += 1
            print('-- running on run# %d (seq.# %d)---------' % (run, runcount))
        
            mastertxt = 'SDSSdata/masterTXT/run%d.txt'%run
            a3d = np.loadtxt(mastertxt, skiprows = 1)
            nfields = int((max(a3d[:, 0])+1) * 6)
            seeing[ii:ii+nfields] = a3d[a3d[:, 2]==args.iBand, 3]
            ii += nfields
        y, bin_edges = np.histogram(seeing, bins=nbins)
        x = (bin_edges[:-1] + bin_edges[1:])/2
        np.savetxt(histtxt, np.transpose(np.vstack((x,y))))
                       
    e = np.sqrt(y)
    e[e<1] = 1
    popt, pcov = optimize.curve_fit(
        lognormal, x, y, sigma=e, p0=[0.5, 0.2, 5e3], absolute_sigma=True)

    # test, SRD values
    # popt[0] = 0.5 #mu
    # popt[1]=0.2 #sigma
    
    yfit = lognormal(x, popt[0], popt[1], popt[2])
    print('mu = %.3f, sigma=%.3f, A=%.1f'%(popt[0], popt[1], popt[2]))
    # plt.errorbar(x,y,e,fmt='ok')
    plt.hist(x,nbins, weights = y)
    plt.plot(x,yfit,'-r')
    plt.grid()
    plt.xlabel('FWHMvk (arcsec)')

    text = 'band = %s'%(bandL[args.iBand])
    plt.annotate(text, xy=(0.5, 0.9),xycoords='axes fraction')
    text = '$\mu$=%.3f $\pm$ %.3f'%(popt[0], np.sqrt(pcov[0,0]))
    plt.annotate(text, xy=(0.5, 0.8),xycoords='axes fraction')
    text = '$\sigma$=%.3f $\pm$ %.3f'%(popt[1], np.sqrt(pcov[1,1]))
    plt.annotate(text, xy=(0.5, 0.7),xycoords='axes fraction')
    # plt.show()
    plt.savefig('output/lognormalHist_%d.png'%args.iBand)

    end = time.time()
    print('time = %8.2fs' % (end - start))
    
if __name__ == "__main__":
    main()
    
