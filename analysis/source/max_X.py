
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize

from sdssinst import sdssinst

"""
find the max and min airmass for each run.
"""


def main():

    objlist = np.loadtxt('data/Stripe82RunList.dat')
    sdss = sdssinst()
    runcount = 0

    nrun = objlist.shape[0]
    runNo = np.zeros(nrun)
    maxX = np.zeros(nrun)
    minX = np.zeros(nrun)
    for line in objlist:

        run = int(line[0])
        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        airmass = txtdata[:, 6]

        runNo[runcount ] = run
        maxX[runcount ] = np.max(airmass)
        minX[runcount ] = np.min(airmass)
        runcount += 1

    plt.figure(figsize=(6,4.5))
    plt.plot(runNo, maxX, 'bo', label='Max Airmass')
    plt.plot(runNo, minX, 'rx', label='Min Airmass')
    plt.legend(loc="upper right", fontsize=11)
    plt.xlabel('Run Number')
    
    plt.tight_layout()
    # plt.show()
    plt.savefig('output/max_X.png')

if __name__ == "__main__":
    main()
