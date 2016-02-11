
# Read a list of (run, minField, maxField) from a text file and
# download psField*fits files in directories
# workingDir/rootName/run/camcol/

# Libraries
import os
import urllib.request
import astropy.io.fits as pyfits
import numpy as np
from matplotlib import pyplot as plt


# 5360 1 99
# psField-005360-1-0099.fit
# http://data.sdss3.org/sas/dr9/boss/photo/redux/301/5360/objcs/1/

# SDSS psf: 2 gaussians plus a power-law core
#  psf = aG1 * exp(-0.5*r^2/sigG1^2) + aG2 * exp(-0.5*r^2/sigG2^2)
#            + aP * (1 + r^2/sigP^2/beta)^(-0.5*beta)
# empirically:
# - aG1 ~ 10*aG2, and sigG2 = ~2*sigG1
# - beta = ~3
# - aP/aG1 = ~10^-3
#  QA plots in above directories (for data, it's easier to use
#                                              photoField files)
# e.g. psPlots1-004874-g[1,2,3].ps

def fetchSDSSphotoField(outdir, run, camcol):
    try:
        infile = "http://data.sdss3.org/sas/dr9/boss/photoObj/301/%d\
/photoField-%06d-%d.fits" % (
            run, run, camcol)
        # Download the fits file
        outfile = outdir + "/%d/photoField-%06d-%d.fits" % (run, run, camcol)
        print("retrieving:%s" % outfile)
        urllib.request.urlretrieve(infile, outfile)
    except:
        print("some problem with run=%06d camCol=%d" % (run, camcol))
        os.remove(outfile)
        pass
    # if the file doesn't exist, urllib still makes an (almost) empty file,
    # remove it...
    statinfo = os.stat(outfile)
    filesize = statinfo.st_size
    if filesize < 300:
        os.remove(outfile)
        return 0
    else:
        print("downloaded photoField-%06d-%d.fits" % (run, camcol))
        return 1


def getSDSSprofRadii():

    radii = np.linspace(0, 15, 16)
    # SDSS profile radii in arcsec, for (more precise) pixel values see p.5 in
    # http://www.astro.princeton.edu/~rhl/photomisc/profiles.ps
    radii[0] = 0.00
    radii[1] = 0.23
    radii[2] = 0.68
    radii[3] = 1.03
    radii[4] = 1.76
    radii[5] = 3.00
    radii[6] = 4.63
    radii[7] = 7.43
    radii[8] = 11.42
    radii[9] = 18.20
    radii[10] = 28.20
    radii[11] = 44.21
    radii[12] = 69.0
    radii[13] = 107.81
    radii[14] = 168.20
    radii[15] = 263.00
    return radii


def plotPSF(data, bandIndex, plotTheory=0):

    # all these are arrays (ugriz), units for sig* are pixels
    sigG1 = data['psf_sigma1'][bandIndex]
    sigG2 = data['psf_sigma2'][bandIndex]
    b = data['psf_b'][bandIndex]
    p0 = data['psf_p0'][bandIndex]
    beta = data['psf_beta'][bandIndex]
    sigP = data['psf_sigmap'][bandIndex]
    nprof = data['prof_nprof'][bandIndex]
    print(sigG1, sigG2, b, p0, beta, sigP, nprof)

    # these are nprof long arrays
    profRadii = getSDSSprofRadii()  # in arcsec
    # better than mean at large radii
    profile = data['prof_med_nmgy'][bandIndex]
    # profile = data['prof_mean_nmgy'][bandIndex]
    profileErr = data['prof_sig_nmgy'][bandIndex]
    # mean root square radius for this annulus (that's how profiles are
    # defined)
    profRadiiMS = np.sqrt(
        (profRadii[:nprof]**2 + profRadii[1:nprof + 1]**2) / 2)
    OKprofRadii = profRadiiMS[:nprof]
    # renormalize to 1 at r~0, and take log10
    # (not exactly at zero because of the mrs radius
    #  but the difference is tiny)
    OKprofile = np.log10(profile[:nprof] / profile[0])
    # error for the log10 of profile
    OKprofileErr = profileErr[:nprof] / profile[0] / np.log(10)

    # best-fit model: double gaussian plus power law
    r = np.linspace(0, 50, 501)
    # for fits, radius must be in pixels
    r2 = (r * 2.5)**2
    psfG1 = np.exp(-0.5 * r2 / (sigG1 * sigG1))
    psfG2 = b * np.exp(-0.5 * r2 / (sigG2 * sigG2))
    # note division by beta! below:
    psfW = p0 * (1 + r2 / (sigP * sigP) / beta)**(-0.5 * beta)
    psfG = psfG1 + psfG2
    # normalized to 1 at r=0 by definition
    psfModel = (psfG + psfW) / (1 + b + p0)
    LpsfG1 = np.log10(psfG1 + 1.0e-9)
    LpsfG2 = np.log10(psfG2)
    LpsfG = np.log10(psfG)
    LpsfW = np.log10(psfW)
    LpsfModel = np.log10(psfModel)

    # get the 2D profile
#    xx,yy = np.meshgrid(r,r)
#    rr2 = (xx*xx+yy*yy) * 2.5**2
#    rrpsfG1 = np.exp(-0.5 * rr2 / (sigG1 * sigG1))
#    rrpsfG2 = b * np.exp(-0.5 * rr2 / (sigG2 * sigG2))
#    # note division by beta! below:
#    rrpsfW = p0 * (1 + rr2 / (sigP * sigP) / beta)**(-0.5 * beta)
#    rrpsfG = rrpsfG1 + rrpsfG2
#    # normalized to 1 at r=0 by definition
#    rrpsfModel = (rrpsfG + rrpsfW) / (1 + b + p0)

    # ------------------------------------------------------------
    # Plot the data and the psf model
    fig1 = plt.figure(figsize=(5, 3.95))
    fig1.subplots_adjust(bottom=0.17, top=0.9, left=0.12,
                         right=0.95, wspace=0.3)
    ax1 = fig1.add_subplot(111)

    ax1.plot(r, LpsfG1, 'k')
    ax1.plot(r, LpsfG2, 'k')
    ax1.plot(r, LpsfG, 'b')
    ax1.plot(r, LpsfW, 'g')
    ax1.plot(r, LpsfModel, 'r', label='PSF fit')

    ax1.errorbar(OKprofRadii, OKprofile, OKprofileErr, fmt='ok')

    ax1.set_xlabel('$r$ (arcsec)')
    ax1.set_ylabel('$PSF$')

    ax1.set_xlim(0, 30.0)
    ax1.set_ylim(-6, 0.5)

    # show which band
    band = {}
    band[0] = "u"
    band[1] = "g"
    band[2] = "r"
    band[3] = "i"
    band[4] = "z"
    text = 'SDSS band: ' + band[bandIndex]
    plt.text(7.0, -0.2, text, fontsize=15, ha='left', va='center')

    # add von Karman predictions from Bo:
    # seeing profile for FWHM=1 arcsec, normalized to 1 at origin
    if (plotTheory):
        # volatile: assumes Bo's file is in local directory!
        Kdata = np.loadtxt('data/r_vonK_Kolm_4.txt',
                           unpack='True')
        radiusK = Kdata[0]
        PK1 = Kdata[1]
        PK2 = Kdata[2]
        PK3 = Kdata[3]
        PK4 = Kdata[4]
        # we will renormalize radius so that the theoretical profile (P1)
        # and best-fit profile agree at the 3 arcsec radius
        # NB theoretical profiles are already normalized to 1 at r=0)
        radiusMatch = 2.0
        # get the value of the best-fit profile at radiusMatch
        profileMatch = np.interp(radiusMatch, r, psfModel)
        # get the radius where theoretical profile is the same as best fit @
        # radiusMatch
        radiusSame1 = np.interp(profileMatch, PK1[::-1], radiusK[::-1])
        # rescale the radius for theoretical profile
        radiusScaled1 = radiusK * (radiusMatch / radiusSame1)
        LPK1 = np.log10(PK1)
        ax1.plot(radiusScaled1, LPK1, 'r--')
        # and now repeat for the second profile
        radiusSame2 = np.interp(profileMatch, PK2[::-1], radiusK[::-1])
        radiusScaled2 = radiusK * (radiusMatch / radiusSame2)
        LPK2 = np.log10(PK2)
        ax1.plot(radiusScaled2, LPK2, 'b-')

        radiusSame3 = np.interp(profileMatch, PK3[::-1], radiusK[::-1])
        radiusScaled3 = radiusK * (radiusMatch / radiusSame3)
        LPK3 = np.log10(PK3)
        ax1.plot(radiusScaled3, LPK3, 'k-')

        radiusSame4 = np.interp(profileMatch, PK4[::-1], radiusK[::-1])
        radiusScaled4 = radiusK * (radiusMatch / radiusSame4)
        LPK4 = np.log10(PK4)
        ax1.plot(radiusScaled4, LPK4, 'm-')

    # figure saved in local directory
    plotName = 'output/SDSSpsf_' + band[bandIndex] + 'Band.png'
    plt.savefig(plotName)

    # plt.show()
    return r, psfModel, OKprofRadii, OKprofile, OKprofileErr


# little helper for plotPSF, e.g.
# # plotPSFfromPhotoField('SDSSdata/4874/photoField-004874-2.fits', 10, 2)
def plotPSFfromPhotoField(filename, ifield, iband):

    hdulist = pyfits.open(filename)
    alldata = hdulist[1].data
    data = alldata[ifield]
    # plot data and best-fit for r band
    radius, psfModel, profRadii, profile, profileErr = plotPSF(data, 2)


# this will download 6 files for 6 camera column from run 4874
if (0):
    outdir = 'SDSSdata'
    run = 4874
    newdir = outdir + "/%d" % (run)
    try:
        os.stat(newdir)
    except:
        os.makedirs(newdir)
    for camcol in range(1, 7):
        print('retrieving file for camcol:', camcol)
        fetchSDSSphotoField(outdir, run, camcol)


# read data and plot PSF for one field from run 4874, camera column 1, and
# all 5 bands
datadir = 'SDSSdata'
run = 4874
camcol = 1
datafile = datadir + "/%d/photoField-%06d-%d.fits" % (run, run, camcol)
hdulist = pyfits.open(datafile)
alldata = hdulist[1].data
for ifield in range(0, 1):
    data = alldata[ifield]
    for iBand in range(1, 2):  # g band only, which includes 500nm
        radius, psfModel, profRadii, profile, profileErr = plotPSF(
            data, iBand, 1)
