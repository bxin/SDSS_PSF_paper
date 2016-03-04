import argparse

import numpy as np
from matplotlib import pyplot as plt

from sdssinst import sdssinst

"""
this code generates one plot for each run.
with run#=-1, it can also project all data from all the runs onto the same plot
but that scatter plot takes quite a long long time to create.
"""


def main():

    parser = argparse.ArgumentParser(
        description='----- fwhm_vk_psfwidth.py ---------')
    parser.add_argument('run', type=int, default=94,
                        help='run number; use -1 for all runs together;\
                        -9 for all runs individually')
    args = parser.parse_args()

    runNo = args.run
    objlist = np.loadtxt('data/Stripe82RunList.dat')
    if runNo > 0:
        # remove all lines but one
        objlist = objlist[objlist[:, 0] == runNo, :]

    sdss = sdssinst()
    runcount = 0

    # vk = np.loadtxt('data/vonK1.0.txt', unpack='True')
    # vk = vk / np.sum(vk)
    # vkneff = 1 / np.sum(vk**2)
    # vkfwhm = 0.664 * 0.1 * np.sqrt(vkneff)
    f, ax1 = plt.subplots(sdss.nBand, sdss.nCamcol, sharex='col',
                          sharey='row', figsize=(12, 8))  # a plot is for a run

    for line in objlist:

        run = int(line[0])
        runcount += 1
        print('-- running on run# %d (seq.# %d)---------' % (run, runcount))

        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        fwhmvk = txtdata[:, 3]
        psf_width = txtdata[:, 4]

        if runNo < -5:
            # a plot is for a run
            f, ax1 = plt.subplots(sdss.nBand, sdss.nCamcol, sharex='col',
                                  sharey='row', figsize=(12, 8))

        for camcol in range(1, sdss.nCamcol + 1):

            for iBand in range(0, sdss.nBand):
                idx = (txtdata[:, 1] == camcol) & (txtdata[:, 2] == iBand)
                ax1[iBand, camcol - 1].scatter(
                    psf_width[idx], fwhmvk[idx], s=2, c='r', edgecolors='r')
                if iBand == 0 and runcount == 1:
                    ax1[iBand, camcol - 1].set_title('camcol=%d' % camcol)
                if camcol == 1 and runcount == 1:
                    text = 'band: %s' % sdss.band[iBand]
                    ax1[iBand, camcol -
                        1].text(1.0, 0.8, text, fontsize=15, ha='left',
                                va='center')
                ax1[iBand, camcol - 1].plot([0.5, 2.5], [0.5, 2.5],
                                            color='k', linestyle='-',
                                            linewidth=2)
                ax1[iBand, camcol - 1].grid()
                if iBand == sdss.nBand - 1:
                    ax1[iBand, camcol - 1].set_xlabel('psf_width')
                if camcol == 1:
                    ax1[iBand, camcol - 1].set_ylabel('FWHMvK')

        if runNo < -5:
            plt.tight_layout()
            plt.savefig(
                'output/fwhmvk_psfwidth/run%d_fwhmvk_psfwidth.png' % (run))
            plt.close()
    if runNo > -5:
        plt.tight_layout()
        if runNo < 0:
            plt.savefig('output/fwhmvk_psfwidth/runAll_fwhmvk_psfwidth.png')
        else:
            plt.savefig(
                'output/fwhmvk_psfwidth/run%d_fwhmvk_psfwidth.png' % (runNo))

        plt.close()

if __name__ == "__main__":
    main()
