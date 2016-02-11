import sys

import numpy as np
from matplotlib import pyplot as plt

from sdssinst import sdssinst


def main():

    objlist = np.loadtxt('data/Stripe82RunList.dat')

    sdss = sdssinst()
    runcount = 0

    vk = np.loadtxt('data/vonK1.0.txt', unpack='True')
    vk = vk / np.sum(vk)
    vkneff = 1 / np.sum(vk**2)
    vkfwhm = 0.664 * 0.1 * np.sqrt(vkneff)
    f, ax1 = plt.subplots(sdss.nBand, sdss.nCamcol, sharex='col',
                          sharey='row', figsize=(12, 8))  # a plot is for a run

    for line in objlist:

        run = int(line[0])
        runcount += 1
        print('-- running on run# %d (seq.# %d)---------' % (run, runcount))

        txtfile = 'SDSSdata/masterTXT/run%d.txt' % (run)
        txtdata = np.loadtxt(txtfile)
        vkscale = txtdata[:, 3]
        psf_width = txtdata[:, 4]

        for camcol in range(1, sdss.nCamcol + 1):

            for iBand in range(0, sdss.nBand):
                ax1[iBand, camcol - 1].plot(psf_width, vkfwhm * vkscale, 'r')
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

        # # plt.suptitle('All units in arcsec ')
        plt.tight_layout()
        # plt.show()
        plt.savefig('output/fwhmvk_psfwidth/run%d_fwhmvk_psfwidth.png' % (run))
        # sys.exit()

    # plt.suptitle('All units in arcsec ')
    plt.tight_layout()
    # plt.show()
    plt.savefig('output/fwhmvk_psfwidth.png')
    sys.exit()

if __name__ == "__main__":
    main()
