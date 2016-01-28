
import os
import time
import argparse

import numpy as np
from matplotlib import pyplot as plt
from sdssinst import sdssinst
from sdssrun import sdssrun
import myTools


def main():
    parser = argparse.ArgumentParser(
        description='----- find_long_wild_runs.py ---------')
    parser.add_argument('vname', help='variable we want to look at, e.g. psf_sigma1,\
psf_sigma2, airmass, neff_psf, psf_nstar, psf_width')
    parser.add_argument('-fitsoff', help='w/o reading in data from fits files',
                        action='store_true')
    parser.add_argument('-tableoff', help='w/o creating the stretchTable',
                        action='store_true')
    args = parser.parse_args()

    start1 = time.time()

    objlist = np.loadtxt('data/Stripe82RunList.dat')

    sdss = sdssinst()

    outTable = 'output/longWildRuns.txt'

    if (not args.tableoff):
        if os.path.isfile(outTable):
            os.remove(outTable)
        fz = open(outTable, 'ab')

        runcount = 0
        for line in objlist:

            start = time.time()
            run = int(line[0])
            runcount += 1
            print('-- running on run# %d (seq.# %d)---------' % (
                run, runcount))

            myRun = sdssrun(run)
            a3dfile = 'output/temp/run%d_%s.txt' % (myRun.runNo, args.vname)
            if (not args.fitsoff):
                a3d = myRun.getBCFtable(sdss, args.vname)
                myTools.savetxt3d(a3dfile, a3d)
            else:
                a3d = np.loadtxt(a3dfile)
                a3d = a3d.reshape(sdss.nBand, sdss.nCamcol, -1)
            longWildTable = myRun.findLongWild(sdss, a3d)
            np.savetxt(fz, longWildTable, fmt='%d %d %d %d %9.6f %9.6f')

            end = time.time()
            print('--- done with run# %d, time elapsed: %8.2fs---' % (
                run, (end - start)))
        fz.close()

    nrun = len(objlist)
    tableData = np.loadtxt(outTable)
    prodbyrun = np.zeros(nrun)
    nfieldsbyrun = np.zeros(nrun)
    pvbyrun = np.zeros(nrun)
    rmsbyrun = np.zeros(nrun)
    runNo = objlist[:, 0].astype(int)
    for i in np.arange(nrun):
        idx = tableData[:, 0] == runNo[i]
        nfieldsbyrun[i] = np.mean(tableData[idx, 1])
        pvbyrun[i] = np.mean(tableData[idx, 4])
        rmsbyrun[i] = np.mean(tableData[idx, 5])
        prodbyrun[i] = nfieldsbyrun[i] * pvbyrun[i] * rmsbyrun[i]
    plt.plot(runNo, nfieldsbyrun,  label='nfields',
             marker='.', color='g', markersize=10)
    plt.plot(runNo, prodbyrun,  label='nfields*PV*RMS (arcsec$^2$)',
             marker='o', color='k', markersize=5)
    plt.plot(runNo, pvbyrun*100,  label='PV (arcsec) x 100',
             marker='v', color='b', markersize=5)
    plt.plot(runNo, rmsbyrun*100,  label='RMS (arcsec) x 100',
             marker='*', color='r', markersize=10)
    plt.xlabel('Run No.')
    plt.ylabel('runs & psf_width variations')
    plt.legend(loc="upper right", fontsize=10)

    plt.savefig('output/longWildRun_%s.png' % (args.vname))
    maxRun = runNo[np.where(prodbyrun == np.max(prodbyrun))]
    runSeq = np.arange(1, nrun + 1)
    for irun in maxRun: #mostly likely there is only 1 number in maxRun
        idx = tableData[:, 0] == irun
        print('longest run with largest %s variation: run %d (seq# %d), \
        nfield = %d, \n\
        averaged over band and camcol: PV = %9.6f, RMS = %9.6f' % (
            args.vname,
            irun,
            runSeq[irun == runNo],
            np.mean(tableData[idx, 1]) ,
            np.mean(tableData[idx, 4]),
            np.mean(tableData[idx, 5])))
    # sys.exit()
    end1 = time.time()
    print('--- Total time elapsed: %8.2fs---' % ((end1 - start1)))

if __name__ == "__main__":
    main()
