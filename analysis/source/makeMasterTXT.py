import time

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

objlist = np.loadtxt('data/Stripe82RunList.dat')

sdss = sdssinst()

start = time.time()

runcount = 0
for line in objlist:

    run = int(line[0])
    runcount += 1
    print('-- running on run# %d (seq.# %d)---------' % (run, runcount))

    myRun = sdssrun(run)
    txtfile = 'SDSSdata/masterTXT/run%d.txt' % (myRun.runNo)
    myRun.writeMasterTXT(sdss, txtfile)
    
end = time.time()
print('time = %8.2fs'%(end-start))
