#!/usr/bin/env python

import os
import numpy as np

def savetxt3d(outfile, a3d):
    if os.path.isfile(outfile):
        os.remove(outfile)
    fz = open(outfile, 'ab')
    for irow in np.arange(0, a3d.shape[0]):
        np.savetxt(fz, a3d[irow,:,:])
        
            
