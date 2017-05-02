import time
import argparse
import multiprocessing

import numpy as np

from sdssrun import sdssrun
from sdssinst import sdssinst

"""
make files, one per run, as a function of field number, with
- field, camera column, bandpass,
- one-parameter fit FWHM from fitting von Karman profile
- other SDSS params: psf_width, airmass, mjd, psf_nstar, neff_psf, sky_frames


E.g. for some run
#  field  camCol  filter   FWHMvK    (other SDSS params)

"""
parser = argparse.ArgumentParser(
    description='-----makeMasterTXT_pool.py------')
parser.add_argument('-p', dest='numproc', default=1, type=int,
                    help='Number of Processors to use')
parser.add_argument('-run', dest='run', type=int, default=-1,
                    help='run number. use -1 to run over all runs. default=-1')
args = parser.parse_args()

objlist = np.loadtxt('data/Stripe82RunList.dat')
# objlist = objlist[0:3,:] #for test

sdss = sdssinst()

start = time.time()


def write1run(argList):

    sdss = argList[0]
    run = argList[1]
    runcount = argList[2]
    print('-- running on run# %d (seq.# %d)---------' % (run, runcount))
    myRun = sdssrun(run)
    txtfile = 'SDSSdata/masterTXT/run%d.txt' % (myRun.runNo)
    myRun.writeMasterTXT(sdss, txtfile)


argList = []
runcount = 0
for line in objlist:

    run = int(line[0])
    if args.run == -1 or args.run == run:
        runcount += 1
        argList.append((sdss, run, runcount))
    # test
    # write1run(argList[0])

pool = multiprocessing.Pool(args.numproc)

pool.map(write1run, argList)
pool.close()
pool.join()

end = time.time()
print('time = %8.2fs' % (end - start))
