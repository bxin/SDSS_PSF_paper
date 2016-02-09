#!/usr/bin/env python

import numpy as np
from scipy import optimize
from scipy import interpolate


class sdsspsf(object):

    def __init__(self, hdu1, ifield, bandIndex):

        data = hdu1[ifield]

        # all these are arrays (ugriz), units for sig* are pixels
        sigG1 = data['psf_sigma1'][bandIndex]
        sigG2 = data['psf_sigma2'][bandIndex]
        b = data['psf_b'][bandIndex]
        p0 = data['psf_p0'][bandIndex]
        beta = data['psf_beta'][bandIndex]
        sigP = data['psf_sigmap'][bandIndex]
        self.nprof = data['prof_nprof'][bandIndex]
        pixScale = data['pixScale'][bandIndex]
        
        # these are nprof long arrays
        self.getSDSSprofRadii()
        # better than mean at large radii
        self.profile = data['prof_med_nmgy'][bandIndex]
        # profile = data['prof_mean_nmgy'][bandIndex]
        self.profileErr = data['prof_sig_nmgy'][bandIndex]
        # mean root square radius for this annulus (that's how profiles are
        # defined)
        profRadiiMS = np.sqrt(
            (self.profRadii[:self.nprof]**2 + self.profRadii[1:self.nprof + 1]**2) / 2)
        self.OKprofRadii = profRadiiMS[:self.nprof]
        # renormalize to 1 at r~0, and take log10
        # (not exactly at zero because of the mrs radius
        #  but the difference is tiny)
        self.OKprofileLinear = self.profile[:self.nprof] / self.profile[0]
        self.OKprofile = np.log10(self.OKprofileLinear)
        # error for the log10 of profile
        self.OKprofileErrLinear = self.profileErr[:self.nprof] / self.profile[0]
        self.OKprofileErr = self.OKprofileErrLinear/ np.log(10)

        # best-fit model: double gaussian plus power law
        self.r = np.linspace(0, 30, 301)
        # for fits, radius must be in pixels
        r2 = (self.r/pixScale)**2
        psfG1 = np.exp(-0.5 * r2 / (sigG1 * sigG1))
        psfG2 = b * np.exp(-0.5 * r2 / (sigG2 * sigG2))
        # note division by beta! below:
        psfW = p0 * (1 + r2 / (sigP * sigP) / beta)**(-0.5 * beta)
        psfG = psfG1 + psfG2
        # normalized to 1 at r=0 by definition
        self.psfModel = (psfG + psfW) / (1 + b + p0)
        LpsfG1 = np.log10(psfG1 + 1.0e-300)  # avoid log10(0)
        LpsfG2 = np.log10(psfG2 + 1.0e-300)
        self.LpsfG = np.log10(psfG + 1.0e-300)
        self.LpsfW = np.log10(psfW + 1.0e-300)
        self.LpsfModel = np.log10(self.psfModel + 1.0e-300)
        #print('p0=%5.3e, sigP=%5.3e, beta=%5.3e' % (p0, sigP, beta))

        # i = self.nprof
        # while (abs(self.OKprofileErrLinear[i-1]/self.OKprofileLinear[i-1]
        #           -0.02)<1e-5 and i-1>=0):
        #     print('i=%d, ratio=%e\n'%(i, self.OKprofileErrLinear[i-1]/self.OKprofileLinear[i-1]))
        #     i -= 1
        # self.nprofErr = i
        self.nprofErr = 4 #i #use 4 points: 0,1,2,3

    def getSDSSprofRadii(self):

        self.profRadii = np.linspace(0, 15, 16)
        # SDSS profile radii in arcsec, for (more precise) pixel values
        # see p.5 in
        # http://www.astro.princeton.edu/~rhl/photomisc/profiles.ps
        self.profRadii[0] = 0.00
        self.profRadii[1] = 0.23
        self.profRadii[2] = 0.68
        self.profRadii[3] = 1.03
        self.profRadii[4] = 1.76
        self.profRadii[5] = 3.00
        self.profRadii[6] = 4.63
        self.profRadii[7] = 7.43
        self.profRadii[8] = 11.42
        self.profRadii[9] = 18.20
        self.profRadii[10] = 28.20
        self.profRadii[11] = 44.21
        self.profRadii[12] = 69.0
        self.profRadii[13] = 107.81
        self.profRadii[14] = 168.20
        self.profRadii[15] = 263.00

    def fit2vonK_curve_fit(self, vonK1arcsec):
        errLinear = self.OKprofileErrLinear.copy()
        errLinear[self.nprofErr:] = 100 
        popt, pcov = optimize.curve_fit(
            lambda r, scaleR, scaleV: scaleVonKR(vonK1arcsec, r, scaleR, scaleV),
            self.OKprofRadii, self.OKprofileLinear, p0=[1.5, 1],
            sigma=errLinear, absolute_sigma=True)

        self.scaleR = popt[0]
        self.scaleV = popt[1]
        print(self.OKprofileErrLinear/self.OKprofileLinear)
        print('scaleR = %7.5f, scaleV=%7.5f\n'%(self.scaleR, self.scaleV))

    def fit2vonK_fminbound(self, vonK1arcsec):
        myfunc = lambda scaleR, scaleV: scaleVonKRChi2(
                vonK1arcsec, self.OKprofRadii, scaleR, scaleV,
            self.OKprofileLinear, self.OKprofileErrLinear)
        xopt = optimize.fminbound(
            myfunc, 0, 3,  xtol=1e-6, maxfun=500, disp=0)

        self.scaleR = xopt[0]
        self.scaleV = xopt[1]


def scaleVonKR(vonK1arcsec, r, scaleR, scaleV):
    vR = scaleR * vonK1arcsec[0, :]
    vv = vonK1arcsec[1, :]
    # stepR = vR[1] - vR[0]
    p = np.zeros(len(r))
    if scaleR > 0:
        f = interpolate.interp1d(vR, vv)
        p = f(r) * scaleV
#        for i in np.arange(len(r)):
#            x1=np.nonzero(vR<r[i])[0][-1]
#            x2=x1+1
#            w1=(vR[x2]-r[i])/stepR
#            w2=(r[i]-vR[x1])/stepR
#            p[i]=vv[x1]*w1+vv[x2]*w2

#    print('scaleR = %7.4f, f0=%7.4f' % (scaleR, p[0]))
    return p


def scaleVonKRChi2(vonK1arcsec, r, scaleR, scaleV, y, err):
    vR = scaleR * vonK1arcsec[0, :]
    vv = vonK1arcsec[1, :]
    stepR = vR[1] - vR[0]
    p = np.zeros(len(r))
    # f = interpolate.interp1d(vR, vv, kind='linear')
    # f = interpolate.interp1d(vR, vv, kind='quadratic')
    f = interpolate.interp1d(vR, vv, kind='cubic')
    p = f(r) * scaleV

    chi2 = np.sum(((p - y)/err)**2)
    # print('---', vonK1arcsec)
    # print(r)
    # print(y)
    # print('scaleR = %7.4f, chi2=%7.4f' % (scaleR, chi2))
    return chi2
