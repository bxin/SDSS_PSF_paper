
import sys
import time
import argparse

import numpy as np
from matplotlib import pyplot as plt
from sdssinst import sdssinst
from sdssrun import sdssrun
import myTools


def main():
    parser = argparse.ArgumentParser(
        description='----- variation_by_run.py ---------')
    parser.add_argument('vname', help='variable we want to look at, e.g. psf_sigma1,\
psf_sigma2, airmass, neff_psf, psf_nstar, psf_width')
    parser.add_argument('irun', type=int, 
                        help='Run Number; use -1 for all runs')
    parser.add_argument('-fitsoff', help='w/o reading in data from fits files',
                        action='store_true')
    args = parser.parse_args()

    start = time.time()
    objlist = np.loadtxt('data/Stripe82RunList.dat')

    sdss = sdssinst()
    totalFields = 0
    runcount = 0
    for line in objlist:
        runcount += 1
        run = int(line[0])
        if ((args.irun>0 ) and (run != args.irun)):
            continue
        print('-- running on run# %d (seq.# %d)---------' % (
            run, runcount))

        f, ax1 = plt.subplots(sdss.nBand, sdss.nCamcol, sharex='col',
                              sharey='row',
                              figsize=(12, 8))  # a plot is for a run
        myRun = sdssrun(run)
        totalFields += myRun.nfields
        myxticks = np.linspace(0, np.ceil(
            myRun.nfields / 100) * 100, np.ceil(myRun.nfields / 100) + 1)
        myxticklabels = ['%d' % (myxticks[i])
                         for i in np.arange(len(myxticks))]
        for i in np.arange(len(myxticks)):
            if i % 2 == 0:
                pass
            else:
                myxticklabels[i] = ''

        if args.vname == 'fwhmvk':
            mastertxt = 'SDSSdata/masterTXT/run%d.txt'%run
            a3d = np.loadtxt(mastertxt, skiprows = 1)
            a3d = a3d[:, 3]
            a3d = a3d.reshape(sdss.nCamcol,-1, sdss.nBand)
            a3d=np.swapaxes(a3d,0,2)
            a3d=np.swapaxes(a3d,1,2)
        else:
            a3dfile = 'output/temp/run%d_%s.txt' % (myRun.runNo, args.vname)
            if (not args.fitsoff):
                a3d = myRun.getBCFtable(sdss, args.vname)
                myTools.savetxt3d(a3dfile, a3d)
            else:
                a3d = np.loadtxt(a3dfile)
                a3d = a3d.reshape(sdss.nBand, sdss.nCamcol, -1)
        ymin = np.min(a3d)
        ymax = np.max(a3d)
            
        for camcol in range(1, sdss.nCamcol + 1):
            for iBand in range(0, sdss.nBand):

                ax1[iBand, camcol - 1].plot(np.arange(0, myRun.nfields),
                                            a3d[iBand, camcol - 1, :], 'b', linewidth=2)
                if iBand == 0:
                    ax1[iBand, camcol - 1].set_title('camcol=%d' % camcol)
                if camcol == 1:
                    text = 'band: %s' % sdss.band[iBand]
                    ax1[iBand, camcol -
                        1].annotate(text, xy=(0.5, 0.8),
                                    xycoords='axes fraction')  # , fontsize=16,
                    # horizontalalignment='right', verticalalignment='bottom')
                text = 'mean=%.3f'%np.mean(a3d[iBand, camcol - 1, :])
                ax1[iBand, camcol -
                        1].annotate(text, xy=(0.2, 0.2),
                                    xycoords='axes fraction')
                if ((args.vname == 'psf_width') or (args.vname == 'fwhmvk')):
                    if iBand==0: #u
                        rat = a3d[iBand, camcol - 1, :-4]/a3d[2, camcol - 1, 4:]
                    elif iBand==1: #g
                        rat = a3d[iBand, camcol - 1, :-8]/a3d[2, camcol - 1, 8:]
                    elif iBand==2: #r
                        rat = a3d[iBand, camcol - 1, :]/a3d[2, camcol - 1, :]
                    elif iBand==3: #g
                        rat = a3d[iBand, camcol - 1, :-2]/a3d[2, camcol - 1, 2:]
                    elif iBand==4: #g
                        rat = a3d[iBand, camcol - 1, :-6]/a3d[2, camcol - 1, 6:]
                    text = 'r=%.2f(%.2f)'%(np.mean(rat), np.std(rat))
                    ax1[iBand, camcol -
                                1].annotate(text, xy=(0.05, 0.05),
                                                xycoords='axes fraction') 
                ax1[iBand, camcol - 1].set_xticks(myxticks)
                ax1[iBand, camcol - 1].set_xticklabels(myxticklabels)
                ax1[iBand, camcol - 1].set_xlim(0, myRun.nfields)
                ax1[iBand, camcol - 1].set_ylim(ymin-(ymax-ymin)*0.5, ymax+(ymax-ymin)*0.5)
                ax1[iBand, camcol - 1].grid()
                
        plt.suptitle('run %d, field %d, %s ' %
                     (run, myRun.nfields, args.vname))
        # plt.tight_layout()
        # plt.show()
        #plt.savefig('output/run%d_%s.png' % (myRun.runNo, args.vname))
        end = time.time()
    print('total number of fields = %d' % totalFields)
    print('time = %8.2fs' % (end - start))
    #    sys.exit()

if __name__ == "__main__":
    main()
