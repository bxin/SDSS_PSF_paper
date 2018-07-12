#!/usr/bin/env python

# import sys
import numpy as np
from scipy import optimize
from scipy import interpolate


class sdsspsf(object):

    def __init__(self, hdu1, ifield, bandIndex, runNo, camcol):

        data = hdu1[ifield]
        self.runNo = runNo
        self.camcol = camcol
        self.field = ifield
        self.band = bandIndex

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
            (self.profRadii[:self.nprof]**2 +
             self.profRadii[1:self.nprof + 1]**2) / 2)
        self.OKprofRadii = profRadiiMS[:self.nprof]
        # renormalize to 1 at r~0, and take log10
        # (not exactly at zero because of the mrs radius
        #  but the difference is tiny)
        self.OKprofileLinear = self.profile[:self.nprof] / self.profile[0]
        self.OKprofile = np.log10(self.OKprofileLinear)
        # error for the log10 of profile
        self.OKprofileErrLinear = self.profileErr[
            :self.nprof] / self.profile[0]
        self.OKprofileErr = self.OKprofileErrLinear / np.log(10)

        # best-fit model: double gaussian plus power law
        self.r = np.linspace(0, 30, 301)
        # for fits, radius must be in pixels
        r2 = (self.r / pixScale)**2
        psfG1 = np.exp(-0.5 * r2 / (sigG1 * sigG1))
        psfG2 = b * np.exp(-0.5 * r2 / (sigG2 * sigG2))
        # note division by beta! below:
        psfW = p0 * (1 + r2 / (sigP * sigP) / beta)**(-0.5 * beta)
        psfG = psfG1 + psfG2
        # normalized to 1 at r=0 by definition
        self.psfModel = (psfG + psfW) / (1 + b + p0)
        # LpsfG1 = np.log10(psfG1 + 1.0e-300)  # avoid log10(0)
        # LpsfG2 = np.log10(psfG2 + 1.0e-300)
        self.LpsfG = np.log10(psfG + 1.0e-300)
        self.LpsfW = np.log10(psfW + 1.0e-300)
        self.LpsfModel = np.log10(self.psfModel + 1.0e-300)
        # print('p0=%5.3e, sigP=%5.3e, beta=%5.3e' % (p0, sigP, beta))

        # i = self.nprof
        # while (abs(self.OKprofileErrLinear[i-1]/self.OKprofileLinear[i-1]
        #           -0.02)<1e-5 and i-1>=0):
        #     print('i=%d, ratio=%e\n'%(i,
        # self.OKprofileErrLinear[i-1]/self.OKprofileLinear[i-1]))
        #     i -= 1
        # self.nprofErr = i
        self.nprofErr = 4  # i #use 4 points: 0,1,2,3

        errLinear = self.OKprofileErrLinear.copy()
        errLinear[self.nprofErr:] = 100
        f = interpolate.interp1d(self.r, self.psfModel, bounds_error=False)
        yy = f(self.OKprofRadii)
        self.G2chi2 = sum(((yy - self.OKprofileLinear) / errLinear)
                          ** 2) / (len(self.OKprofRadii) - 2)
        idx = self.OKprofRadii < 2
        self.G2chi2lr = sum(
            ((yy[idx] - self.OKprofileLinear[idx]) / errLinear[idx])**2) / (sum(idx))
        idx = self.OKprofRadii >= 2
        self.G2chi2hr = sum(
            ((yy[idx] - self.OKprofileLinear[idx]) / errLinear[idx])**2) / sum(idx)

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

    def fit2vonK_curve_fit(self, vonK1arcsec, vonK2D=None, grid1d=None):
        """
        vonk2D is only used for verification when we do the convolution.
        The fit for scaleR and scaleV is done in 1D.
        """

        errLinear = self.OKprofileErrLinear.copy()
        errLinear[self.nprofErr:] = 100
        try:
            popt, pcov = optimize.curve_fit(
                lambda r, scaleR, scaleV: scaleVonKR(
                    vonK1arcsec, r, scaleR, scaleV),
                self.OKprofRadii, self.OKprofileLinear, p0=[1.5, 1],
                sigma=errLinear, absolute_sigma=True)

            self.scaleR = popt[0]
            self.scaleV = popt[1]
            # for curve with von Karman only
            self.vvR = vonK1arcsec[0] * self.scaleR
            self.vvv = vonK1arcsec[1] * self.scaleV            
            if not (vonK2D is None):
                newN = convVonK2D(vonK2D, grid1d, self.scaleR,
                                  self.scaleV, 1, self.tailP)
                m = int((len(grid1d) - 1) / 2)
                self.vR = grid1d[m:]
                m1 = max(np.argmax(newN == np.max(newN), axis=0))
                m2 = max(np.argmax(newN == np.max(newN), axis=1))
                self.vv = newN[m1, m2:m2 + m + 1]
            else:
                self.vR = self.vvR
                self.vv = self.vvv
            yy = scaleVonKR(vonK1arcsec, self.OKprofRadii,
                            self.scaleR, self.scaleV)
            self.chi2 = sum(((yy - self.OKprofileLinear) / errLinear)
                            ** 2) / (len(self.OKprofRadii) - 2)
            idx = self.OKprofRadii < 2
            self.chi2lr = sum(
                ((yy[idx] - self.OKprofileLinear[idx]) / errLinear[idx])**2) / (sum(idx))
            idx = self.OKprofRadii >= 2
            self.chi2hr = sum(
                ((yy[idx] - self.OKprofileLinear[idx]) / errLinear[idx])**2) / sum(idx)

        except (RuntimeError, ValueError) as e:
            print('in fit2vonK_curve_fit\n')
            print(e)
            print('run#=%d, camcol=%d, field=%d, band=%d\n' % (
                self.runNo, self.camcol, self.field, self.band))
            print('x=\n')
            print(self.OKprofRadii)
            print('y=\n')
            print(self.OKprofileLinear)
            print('err=\n')
            print(errLinear)
            self.scaleR = -999
            self.scaleV = -999
            # sys.exit()
            
    def fit2conv_curve_fit(self, vonK2D, grid1d):
        # use gaussian+10^[eta*(ax^2+bx+c)] as instrumental PSF 
        try:
            errLinear = self.OKprofileErrLinear.copy()
            errLinear[self.nprofErr:] = 100
            assert self.scaleR>-1, "Exiting fit2conv_curve_fit, scaleR<-1" 
            idx = min(abs(self.OKprofRadii - 15)) == abs(self.OKprofRadii - 15)
            errLinear[idx] = 1e-10
            popt, pcov = optimize.curve_fit(
                lambda r, eta: convVonK(
                    vonK2D, grid1d, r, self.scaleR, self.scaleV, eta, self.tailP),
                self.OKprofRadii, self.OKprofileLinear, p0=[1],
                bounds = (0.1, 10), 
                sigma=errLinear, absolute_sigma=True)

            self.tailEta = popt[0]

            newN = convVonK2D(vonK2D, grid1d, self.scaleR,
                              self.scaleV, self.tailEta, self.tailP)
            m = int((len(grid1d) - 1) / 2)
            self.vR = grid1d[m:]
            m1 = max(np.argmax(newN == np.max(newN), axis=0))
            m2 = max(np.argmax(newN == np.max(newN), axis=1))
            self.vv = newN[m1, m2:m2 + m + 1]

            # get the tail
            sigma = self.tailP[2]
            x1= self.tailP[3]
            y1= self.tailP[4]
            x2 = self.tailP[5]
            y2 = self.tailP[6]
            x3 = self.tailP[7]
            y3 = self.tailP[8]
            abcM = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
            abc = np.dot(np.linalg.inv(abcM), np.array([[y1] ,[y2], [y3]]))
            G1 = 10**(self.tailEta*(abc[0][0]*self.vR*self.vR+abc[1][0]*self.vR+abc[2][0]))

            self.vtail = np.exp(-(self.vR**2) / 2 / sigma**2) + G1
    
            yy = convVonK(vonK2D, grid1d, self.OKprofRadii, self.scaleR, self.scaleV,
                          self.tailEta, self.tailP)
            self.chi2 = sum(((yy - self.OKprofileLinear) / errLinear)
                            ** 2) / (len(self.OKprofRadii) - 2)
            idx = self.OKprofRadii < 2
            self.chi2lr = sum(
                ((yy[idx] - self.OKprofileLinear[idx]) / errLinear[idx])**2) / (sum(idx))
            idx = self.OKprofRadii >= 2
            self.chi2hr = sum(
                ((yy[idx] - self.OKprofileLinear[idx]) / errLinear[idx])**2) / sum(idx)

        except (RuntimeError, ValueError, AssertionError) as e:
            print('in fit2conv_curve_fit\n')
            print(e)
            print('run#=%d, camcol=%d, field=%d, band=%d\n' % (
                self.runNo, self.camcol, self.field, self.band))
            print('x=\n')
            print(self.OKprofRadii)
            print('y=\n')
            print(self.OKprofileLinear)
            print('err=\n')
            print(errLinear)
            self.tailEta= -999
            if self.scaleR < -1:
                self.chi2 = -999
                self.chi2lr = -999
                self.chi2hr = -999
                self.vvR = -999
                self.vvv = -999
                self.vR = -999
                self.vv = -999
            
            # sys.exit()
            
    def fit2conv_curve_fit_log(self, vonK2D, grid1d):
        # use gaussian+10^[eta*(ax^2+bx+c)] as instrumental PSF 
        try:
            errLog = self.OKprofileErr.copy()
            # errLinear[self.nprofErr:] = 100
            assert self.scaleR>-1, "Exiting fit2conv_curve_fit_log, scaleR<-1" 
            # idx = min(abs(self.OKprofRadii - 15)) == abs(self.OKprofRadii - 15)
            # errLinear[idx] = 1e-10
            popt, pcov = optimize.curve_fit(
                lambda r, eta: logConvVonK(
                    vonK2D, grid1d, r, self.scaleR, self.scaleV, eta, self.tailP),
                self.OKprofRadii, self.OKprofile, p0=[1],
                bounds = (0.1, 10), 
                sigma=errLog, absolute_sigma=True)

            self.tailEta = popt[0]

            newN = convVonK2D(vonK2D, grid1d, self.scaleR,
                              self.scaleV, self.tailEta, self.tailP)
            m = int((len(grid1d) - 1) / 2)
            self.vR = grid1d[m:]
            m1 = max(np.argmax(newN == np.max(newN), axis=0))
            m2 = max(np.argmax(newN == np.max(newN), axis=1))
            self.vv = newN[m1, m2:m2 + m + 1]

            yy = logConvVonK(vonK2D, grid1d, self.OKprofRadii, self.scaleR, self.scaleV,
                          self.tailEta, self.tailP)
            self.chi2 = sum(((yy - self.OKprofile) / errLog)
                            ** 2) / (len(self.OKprofRadii) - 2)
            idx = self.OKprofRadii < 2
            self.chi2lr = sum(
                ((yy[idx] - self.OKprofile[idx]) / errLog[idx])**2) / (sum(idx))
            idx = self.OKprofRadii >= 2
            self.chi2hr = sum(
                ((yy[idx] - self.OKprofile[idx]) / errLog[idx])**2) / sum(idx)

        except (RuntimeError, ValueError, AssertionError) as e:
            print('in fit2conv_curve_fit\n')
            print(e)
            print('run#=%d, camcol=%d, field=%d, band=%d\n' % (
                self.runNo, self.camcol, self.field, self.band))
            print('x=\n')
            print(self.OKprofRadii)
            print('y=\n')
            print(self.OKprofileLinear)
            print('err=\n')
            print(errLinear)
            self.tailEta= -999
            if self.scaleR < -1:
                self.chi2 = -999
                self.chi2lr = -999
                self.chi2hr = -999
                self.vvR = -999
                self.vvv = -999
                self.vR = -999
                self.vv = -999
            
            # sys.exit()

    def fit2conv_curve_fit_log_ab(self, vonK2D, grid1d):
        # use gaussian+10^[eta*(ax^2+bx+c)] as instrumental PSF 
        try:
            errLog = self.OKprofileErr.copy()
            # errLinear[self.nprofErr:] = 100
            assert self.scaleR>-1, "Exiting fit2conv_curve_fit_log, scaleR<-1" 
            # idx = min(abs(self.OKprofRadii - 15)) == abs(self.OKprofRadii - 15)
            # errLinear[idx] = 1e-10
            popt, pcov = optimize.curve_fit(
                lambda r, eta: logConvVonKAB(
                    vonK2D, grid1d, r, self.scaleR, self.scaleV, eta, self.tailP[2], self.tailP[3]),
                self.OKprofRadii, self.OKprofile, p0=[self.tailP[4]],
                bounds = (min(self.tailP[4]*10, self.tailP[4]*0.1), max(self.tailP[4]*10, self.tailP[4]*0.1)), 
                sigma=errLog, absolute_sigma=True)

            self.tailEta = popt[0]

            newN = convVonKAB2D(vonK2D, grid1d, self.scaleR,
                              self.scaleV, self.tailEta, self.tailP[2], self.tailP[3])
            m = int((len(grid1d) - 1) / 2)
            self.vR = grid1d[m:]
            m1 = max(np.argmax(newN == np.max(newN[100:-100,100:-100]), axis=0))
            m2 = max(np.argmax(newN == np.max(newN[100:-100,100:-100]), axis=1))
            self.vv = newN[m1, m2:m2 + m + 1]

            G1 = 10**(self.tailEta*(self.tailP[2]*self.vR*self.vR+self.tailP[3]*self.vR+1))
            sigma = 0.1
            self.vtail = np.exp(-(self.vR**2) / 2 / sigma**2) + G1
            
            yy = logConvVonKAB(vonK2D, grid1d, self.OKprofRadii, self.scaleR, self.scaleV,
                          self.tailEta, self.tailP[2], self.tailP[3])
            self.chi2 = sum(((yy - self.OKprofile) / errLog)
                            ** 2) / (len(self.OKprofRadii) - 2)
            idx = self.OKprofRadii < 2
            self.chi2lr = sum(
                ((yy[idx] - self.OKprofile[idx]) / errLog[idx])**2) / (sum(idx))
            idx = self.OKprofRadii >= 2
            self.chi2hr = sum(
                ((yy[idx] - self.OKprofile[idx]) / errLog[idx])**2) / sum(idx)

        except (RuntimeError, ValueError, AssertionError) as e:
            print('in fit2conv_curve_fit\n')
            print(e)
            print('run#=%d, camcol=%d, field=%d, band=%d\n' % (
                self.runNo, self.camcol, self.field, self.band))
            print('x=\n')
            print(self.OKprofRadii)
            print('y=\n')
            print(self.OKprofileLinear)
            print('err=\n')
            print(errLinear)
            self.tailEta= -999
            if self.scaleR < -1:
                self.chi2 = -999
                self.chi2lr = -999
                self.chi2hr = -999
                self.vvR = -999
                self.vvv = -999
                self.vR = -999
                self.vv = -999
            
            # sys.exit()
            
    def fit2conv_curve_fit_fitab(self, vonK2D, grid1d):
        # use gaussian+10^[eta*(ax^2+bx+c)] as instrumental PSF 
        try:
            errLog = self.OKprofileErr.copy()
            # errLinear[self.nprofErr:] = 100
            assert self.scaleR>-1, "Exiting fit2conv_curve_fit_ab, scaleR<-1" 
            # idx = min(abs(self.OKprofRadii - 15)) == abs(self.OKprofRadii - 15)
            # errLinear[idx] = 1e-10
            x1= self.tailP[3]
            y1= self.tailP[4]
            x2 =self.tailP[5]
            y2 =self.tailP[6]
            x3 =self.tailP[7]
            y3 =self.tailP[8]
            abcM = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
            abc = np.dot(np.linalg.inv(abcM), np.array([[y1] ,[y2], [y3]]))
            tailA0 = abc[0][0]/abc[2][0]
            tailB0 = abc[1][0]/abc[2][0]
            tailEta0 = abc[2][0]
            popt, pcov = optimize.curve_fit(
                lambda r, eta, tailA, tailB: logConvVonKAB(
                    vonK2D, grid1d, r, self.scaleR, self.scaleV, eta, tailA, tailB),
                self.OKprofRadii, self.OKprofile, p0=[tailEta0, tailA0, tailB0],
                bounds = ([min(tailEta0*10, tailEta0*0.1),-2*abs(tailA0),-2*abs(tailB0)],[max(tailEta0*10, tailEta0*0.1),2*abs(tailA0),2*abs(tailB0)]),
                sigma=errLog, absolute_sigma=True)

            self.tailEta = popt[0]
            self.tailA = popt[1]
            self.tailB = popt[2]
            print('From fit: tailA = %.1e, tailB = %.1e\n'%(self.tailA, self.tailB))
            print('uncertainties: tailA = %.1e, tailB = %.1e\n'%(np.sqrt(pcov[1,1]), np.sqrt(pcov[2,2])))
            
            newN = convVonKAB2D(vonK2D, grid1d, self.scaleR,
                              self.scaleV, self.tailEta, self.tailA, self.tailB)
            m = int((len(grid1d) - 1) / 2)
            self.vR = grid1d[m:]
            m1 = max(np.argmax(newN == np.max(newN), axis=0))
            m2 = max(np.argmax(newN == np.max(newN), axis=1))
            self.vv = newN[m1, m2:m2 + m + 1]

            G12d = 10**(self.tailEta*(self.tailA*self.vR*self.vR+self.tailB*self.vR+1))
            sigma = 0.1
            self.vtail = np.exp(-(self.vR**2) / 2 / sigma**2) + G12d
            
            yy = logConvVonKAB(vonK2D, grid1d, self.OKprofRadii, self.scaleR, self.scaleV,
                          self.tailEta, self.tailA, self.tailB)
            self.chi2 = sum(((yy - self.OKprofile) / errLog)
                            ** 2) / (len(self.OKprofRadii) - 2)
            idx = self.OKprofRadii < 2
            self.chi2lr = sum(
                ((yy[idx] - self.OKprofile[idx]) / errLog[idx])**2) / (sum(idx))
            idx = self.OKprofRadii >= 2
            self.chi2hr = sum(
                ((yy[idx] - self.OKprofile[idx]) / errLog[idx])**2) / sum(idx)

        except (RuntimeError, ValueError, AssertionError) as e:
            print('in fit2conv_curve_fit\n')
            print(e)
            print('run#=%d, camcol=%d, field=%d, band=%d\n' % (
                self.runNo, self.camcol, self.field, self.band))
            print('x=\n')
            print(self.OKprofRadii)
            print('y=\n')
            print(self.OKprofileLinear)
            print('err=\n')
            print(errLog)
            self.tailEta= -999
            if self.scaleR < -1:
                self.chi2 = -999
                self.chi2lr = -999
                self.chi2hr = -999
                self.vvR = -999
                self.vvv = -999
                self.vR = -999
                self.vv = -999
            
            # sys.exit()
            
    # used to use fminbound(), but it is for 1D optimization
    # tried fmin_cg(), but our gradient can only be calculated numerically
    # fmin() takes longer and should be more accurate,
    # we only use it when curve_fit() fails.
    def fit2vonK_fmin(self, vonK1arcsec, vonK2D=None, grid1d=None):
        """
        vonk2D is only used for verification when we do the convolution.
        The fit for scaleR and scaleV is done in 1D.
        """

        errLinear = self.OKprofileErrLinear.copy()
        errLinear[self.nprofErr:] = 100

        try:
            # raise RuntimeError
            xopt = optimize.fmin(
                lambda scaleRV: scaleVonKRChi2(
                    vonK1arcsec, self.OKprofRadii, scaleRV,
                    self.OKprofileLinear, errLinear),
                [1.5, 1], disp=1)

            self.scaleR = xopt[0]
            self.scaleV = xopt[1]
            self.vvR = vonK1arcsec[0] * self.scaleR
            self.vvv = vonK1arcsec[1] * self.scaleV
            
            if not (vonK2D is None):
                newN = convVonK2D(vonK2D, grid1d, self.scaleR,
                                  self.scaleV, 1, self.tailP)
                m = int((len(grid1d) - 1) / 2)
                self.vR = grid1d[m:]
                m1 = max(np.argmax(newN == np.max(newN), axis=0))
                m2 = max(np.argmax(newN == np.max(newN), axis=1))
                self.vv = newN[m1, m2:m2 + m + 1]
            else:
                self.vR = self.vvR
                self.vv = self.vvv

            self.chi2 = (scaleVonKRChi2(
                vonK1arcsec, self.OKprofRadii, xopt, self.OKprofileLinear, errLinear)) \
                / (len(self.OKprofRadii) - 2)
            idx = self.OKprofRadii < 2
            self.chi2lr = (scaleVonKRChi2(
                vonK1arcsec, self.OKprofRadii[idx], xopt, self.OKprofileLinear[idx], errLinear[idx])) / (sum(idx))
            idx = self.OKprofRadii >= 2
            self.chi2hr = (scaleVonKRChi2(
                vonK1arcsec, self.OKprofRadii[idx], xopt, self.OKprofileLinear[idx], errLinear[idx])) / sum(idx)

        except (RuntimeError, ValueError) as e:
            print('in fit2vonK_fmin\n')
            print(e)
            print('run#=%d, camcol=%d, field=%d, band=%d\n' % (
                self.runNo, self.camcol, self.field, self.band))
            print('x=\n')
            print(self.OKprofRadii)
            print('y=\n')
            print(self.OKprofileLinear)
            print('err=\n')
            print(errLinear)
            self.scaleR = -999
            self.scaleV = -999
            # sys.exit()


def scaleVonKR(vonK1arcsec, r, scaleR, scaleV):

    vR = scaleR * vonK1arcsec[0, :]
    vv = vonK1arcsec[1, :]
    # stepR = vR[1] - vR[0]
    p = np.zeros(len(r))
    if scaleR > 0:
        f = interpolate.interp1d(vR, vv, bounds_error=False)
        p = f(r) * scaleV
#        for i in np.arange(len(r)):
#            x1=np.nonzero(vR<r[i])[0][-1]
#            x2=x1+1
#            w1=(vR[x2]-r[i])/stepR
#            w2=(r[i]-vR[x1])/stepR
#            p[i]=vv[x1]*w1+vv[x2]*w2

#    print('scaleR = %7.4f, f0=%7.4f' % (scaleR, p[0]))
    return p


def scaleVonKRChi2(vonK1arcsec, r, scaleRV, y, err):

    scaleR = scaleRV[0]
    scaleV = scaleRV[1]
    vR = scaleR * vonK1arcsec[0, :]
    vv = vonK1arcsec[1, :]
    # stepR = vR[1] - vR[0]
    p = np.zeros(len(r))
    # f = interpolate.interp1d(vR, vv, kind='linear')
    # f = interpolate.interp1d(vR, vv, kind='quadratic')
    f = interpolate.interp1d(vR, vv, kind='cubic', bounds_error=False)
    p = f(r) * scaleV

    chi2 = np.sum(((p - y) / err)**2)
    # print('---', vonK1arcsec)
    # print(r)
    # print(y)
    # print('scaleR = %7.4f, chi2=%7.4f' % (scaleR, chi2))
    return chi2



def convVonK(vonK2D, grid1d, r, scaleR, scaleV, eta, tailP):

    newN = convVonK2D(vonK2D, grid1d, scaleR, scaleV, eta, tailP)
    m = int((len(grid1d) - 1) / 2)
    vR = grid1d[m:]
    m1 = max(np.argmax(newN == np.max(newN), axis=0))
    m2 = max(np.argmax(newN == np.max(newN), axis=1))
    vv = newN[m1, m2:m2 + m + 1]

    # stepR = vR[1] - vR[0]
    p = np.zeros(len(r))
    if scaleR > 0:
        f = interpolate.interp1d(vR, vv, bounds_error=False)
        p = f(r)

        # print('scaleR = %7.4f, eta=%7.4f, p(15arcsec)=%e' % (
        #     scaleR, eta, f(15)))
    return p

def convVonKAB(vonK2D, grid1d, r, scaleR, scaleV, eta, tailA, tailB):

    newN = convVonKAB2D(vonK2D, grid1d, scaleR, scaleV, eta, tailA, tailB)
    m = int((len(grid1d) - 1) / 2)
    vR = grid1d[m:]
    m1 = max(np.argmax(newN == np.max(newN[100:-100,100:-100]), axis=0))
    m2 = max(np.argmax(newN == np.max(newN[100:-100,100:-100]), axis=1))
    vv = newN[m1, m2:m2 + m + 1]
    # print('m1=%d, m2=%d, eta=%e, tailA=%e, tailB=%e\n'%(m1, m2, eta, tailA, tailB))
    
    # stepR = vR[1] - vR[0]
    p = np.zeros(len(r))
    if (scaleR > 0) and (vv.shape[0]==vR.shape[0]):
        f = interpolate.interp1d(vR, vv, bounds_error=False)
        p = f(r)

        # print('scaleR = %7.4f, eta=%7.4f, p(15arcsec)=%e' % (
        #     scaleR, eta, f(15)))
    # print(np.any(~np.isfinite(vv)))
    return p

def logConvVonK(vonK2D, grid1d, r, scaleR, scaleV, eta, tailP):
    return np.log10(convVonK(vonK2D, grid1d, r, scaleR, scaleV, eta, tailP))

def logConvVonKAB(vonK2D, grid1d, r, scaleR, scaleV, eta, tailA, tailB):
    aaa = convVonKAB(vonK2D, grid1d, r, scaleR, scaleV, eta, tailA, tailB)
    # print(aaa)
    if np.any(aaa==0):
        return aaa
    else:
        return np.log10(aaa)

def convVonK2D(vonK2D, grid1d, scaleR, scaleV, eta, tailP):

    x1 = grid1d * scaleR
    f = interpolate.RectBivariateSpline(x1, x1, vonK2D, kx=1, ky=1)
    vk = f(grid1d, grid1d)
    vkfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(vk), s=vk.shape))

    xm, ym = np.meshgrid(grid1d, grid1d)
    rm = np.sqrt(xm*xm+ym*ym)
    sigma = tailP[2]
    x1= tailP[3]
    y1= tailP[4]
    x2 = tailP[5]
    y2 = tailP[6]
    x3 = tailP[7]
    y3 = tailP[8]
    abcM = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
    abc = np.dot(np.linalg.inv(abcM), np.array([[y1] ,[y2], [y3]]))
    G12d = 10**(eta*(abc[0][0]*rm*rm+abc[1][0]*rm+abc[2][0]))

    g2d = np.exp(-(xm * xm + ym * ym) / 2 / sigma**2) + G12d
    psffft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g2d), s=g2d.shape))
    prodfft = psffft * vkfft

    new = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                   s=prodfft.shape)))
    newN = new / np.max(new) * scaleV
    return newN

def convVonKAB2D(vonK2D, grid1d, scaleR, scaleV, eta, tailA, tailB):

    x1 = grid1d * scaleR
    f = interpolate.RectBivariateSpline(x1, x1, vonK2D, kx=1, ky=1)
    vk = f(grid1d, grid1d)
    vkfft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(vk), s=vk.shape))

    xm, ym = np.meshgrid(grid1d, grid1d)
    rm = np.sqrt(xm*xm+ym*ym)
    G12d = 10**(eta*(tailA*rm*rm+tailB*rm+1))

    sigma = 0.1
    g2d = np.exp(-(xm * xm + ym * ym) / 2 / sigma**2) + G12d
    psffft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(g2d), s=g2d.shape))
    prodfft = psffft * vkfft

    new = np.absolute(np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(prodfft),
                                                   s=prodfft.shape)))
    newN = new / np.max(new) * scaleV
    return newN
