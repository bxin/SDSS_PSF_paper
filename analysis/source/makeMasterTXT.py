import time
import multiprocessing
import argparse

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
    description='-----makeMasterTXT.py------')
parser.add_argument('-p', dest='numproc', default=1, type=int,
                    help='Number of Processors Phosim uses')
args = parser.parse_args()

objlist = np.loadtxt('data/Stripe82RunList.dat')
objlist = objlist[0:2,:]

sdss = sdssinst()

start = time.time()

def write1run(sdss, run, runcount):
    print('-- running on run# %d (seq.# %d)---------' % (run, runcount))
    myRun = sdssrun(run)
    txtfile = 'SDSSdata/masterTXT/run%d.txt' % (myRun.runNo)
    myRun.writeMasterTXT(sdss, txtfile)
    

jobs = []
counter = 0

runcount = 0
for line in objlist:

    run = int(line[0])
    runcount += 1
    p = multiprocessing.Process(
        target=write1run, args=(sdss, run, runcount))
    jobs.append(p)
    p.start()
    counter += 1
    if (counter == args.numproc) or (run == objlist[-1, 0]):
        for p in jobs:
            p.join()
        counter = 0
        jobs = []
    
end = time.time()
print('time = %8.2fs'%(end-start))

