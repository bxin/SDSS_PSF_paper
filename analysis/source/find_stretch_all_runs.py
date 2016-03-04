
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
        description='----- find_stretch_all_runs.py ---------')
    parser.add_argument('vname', help='variable we want to look at, e.g. psf_sigma1,\
psf_sigma2, airmass, neff_psf, psf_nstar, psf_width')
    parser.add_argument('-pvfrac', type=float, default=0.1,
                        help='tolerance on peak-to-peak variation, \
default=0.1')
    parser.add_argument('-rmsfrac', type=float, default=0.05,
                        help='tolerance on RMS of the stretch, default=0.05')
    parser.add_argument('-fitsoff', help='w/o reading in data from fits files',
                        action='store_true')
    parser.add_argument('-tableoff', help='w/o creating the stretchTable',
                        action='store_true')
    args = parser.parse_args()

    start1 = time.time()

    objlist = np.loadtxt('data/Stripe82RunList.dat')

    sdss = sdssinst()

    outTable = 'output/stretchTable.txt'

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
            stretchTable = myRun.findStretch(
                sdss, a3d, args.pvfrac, args.rmsfrac)
            np.savetxt(fz, stretchTable, fmt='%d')

            end = time.time()
            print('--- done with run# %d, time elapsed: %8.2fs---' % (
                run, (end - start)))
        fz.close()

    nrun = len(objlist)
    tableData = np.loadtxt(outTable)
    maxbyrun = np.zeros(nrun)
    runNo = objlist[:, 0].astype(int)
    for i in np.arange(nrun):
        # last column in tableTable is the length of the stretch
        maxbyrun[i] = np.max(tableData[tableData[:, 0] == runNo[i], -1])
    plt.plot(runNo, maxbyrun,  # label='',
             marker='o', color='r', markersize=10)
    plt.xlabel('Run No.')
    plt.ylabel('Longest Stretch')

    plt.savefig('output/maxStretch_%s.png' % (args.vname))
    for iRecord in range(0,2): #number of records to print out
        maxRow = np.where(tableData[:, -1] == np.max(tableData[:, -1]))
        runSeq = np.arange(1, nrun + 1)
        for irow in maxRow[0]:
            irun = tableData[irow, 0]
            print('longest stretch found in run %d (seq# %d), \
    field range = (%d, %d), %d runs' % (
                irun,
                runSeq[irun == runNo],
                tableData[irow, 1], tableData[irow, 2], tableData[irow, 3]))
        rmRow = np.where(tableData[:, 0] == irun)
        tableData[rmRow[0], -1] = 0
    # sys.exit()
    end1 = time.time()
    print('--- Total time elapsed: %8.2fs---' % ((end1 - start1)))

if __name__ == "__main__":
    main()
