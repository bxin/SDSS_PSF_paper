#!/usr/bin/env python

import numpy as np

from astropy.io import fits
from sdsspsf import sdsspsf

class sdssrun(object):

    def __init__(self, runNo):
        rootName = "SDSSdata"
        self.runNo = runNo
        self.fitsdir = "%s/%d/" % (rootName, runNo)
        
        datafile = self.fitsdir + "photoField-%06d-%d.fits" % (self.runNo, 1)
        hdulist = fits.open(datafile)
        self.nfields = hdulist[0].header['NFIELDS']
            
    def getBCFtable(self, sdss, vname):
        if vname == 'vkfwhm':
            vdata = np.loadtxt('data/r_vonK_Kolm.txt',
                                unpack='True')
            radius = vdata[0]
            vonK = vdata[1]
            vonK1arcsec = np.vstack((radius, vonK))
        
        a3d = np.zeros((sdss.nBand, sdss.nCamcol, self.nfields))
        for camcol in range(1, sdss.nCamcol + 1):
            print('running on camcol#%d' % camcol)
            datafile = self.fitsdir + "photoField-%06d-%d.fits" % (self.runNo, camcol)

            hdulist = fits.open(datafile)
            hdu1 = hdulist[1].data
            for ifield in range(0, self.nfields):
                if ifield % 100 == 0:
                    print('field No. = %d/%d' % (ifield, self.nfields))
                for iBand in range(0, sdss.nBand):
                    if vname == 'vkfwhm':
                        psf = sdsspsf(hdu1, ifield, iBand)
                        psf.fit2vonK_curve_fit(vonK1arcsec)
                        a3d[iBand, camcol - 1, ifield] = psf.scaleR
                    else:
                        # a3d[iBand, camcol - 1, ifield] =
                        # getattr(psf, args.vname)
                        a3d[iBand, camcol - 1,
                            ifield] = hdu1[ifield][vname][iBand]
        return a3d

    def findStretch(self, sdss, a3d, pvfrac, rmsfrac):
        dataTable = np.zeros((self.nfields, 4))
        for ifield in range(0, self.nfields):
            endfield = ifield
            for jfield in range(ifield + 1, self.nfields):
                goodStretch = True
                for camcol in range(1, sdss.nCamcol + 1):
                    for iBand in range(0, sdss.nBand):
                        stretch = a3d[iBand, camcol - 1, ifield:jfield + 1]
                        ave = np.mean(stretch)
                        # if (np.max(stretch)-ave)/ave > pvfrac:
                        #     goodStretch = False
                        # if (ave - np.min(stretch))/ave > pvfrac:
                        #     goodStretch = False
                        if (np.max(stretch) - np.min(stretch)) / ave \
                                > pvfrac:
                            goodStretch = False
                            break
                        if np.std(stretch) / ave > rmsfrac:
                            goodStretch = False
                            break
                    if (not goodStretch):
                        break
                if goodStretch:
                    endfield = jfield
                else:
                    break
            dataTable[ifield, :] = [self.runNo, ifield, endfield, endfield - ifield + 1]
        return dataTable

    def findLongWild(self, sdss, a3d):
        # 6 columns: run#, nfields, band, camcol, pv, rms
        dataTable = np.zeros((sdss.nBand*sdss.nCamcol, 6))
        pv = np.ptp(a3d, axis=2) #5 x 6 np array
        rms = np.std(a3d, axis=2) #5 x 6 np array
        for camcol in range(1, sdss.nCamcol + 1):
            for iBand in range(0, sdss.nBand):
                iRow = iBand*sdss.nCamcol + camcol - 1
                dataTable[iRow, :] = [self.runNo, self.nfields, iBand, camcol,
                                     pv[iBand , camcol-1], rms[iBand, camcol-1] ]
        
        return dataTable
        
    def writeMasterTXT(self, sdss, filename):

        vdata = np.loadtxt('data/r_vonK_Kolm.txt',
                            unpack='True')
        radius = vdata[0]
        vonK = vdata[1]
        vonK1arcsec = np.vstack((radius, vonK))
        
        fid = open(filename, 'w')
        fid.write('#field \t camCol  filter FWHMvK ')
        otherParams=[['psf_width', '%5.3f'],
                     ['airmass', '%5.3f'],
                     ['mjd\t\t', '%12.6f'],
                     ['psf_nstar ', '%d'],
                     ['neff_psf\t', '%7.3f'],
                     ['sky_frames\t', '%5.3f']]
        for i in range(len(otherParams)):
            fid.write('%s '%otherParams[i][0])
        fid.write('\n')
        
        for camcol in range(1, sdss.nCamcol + 1):
            print('running on camcol#%d' % camcol)
            datafile = self.fitsdir + "photoField-%06d-%d.fits" % (self.runNo, camcol)

            hdulist = fits.open(datafile)
            hdu1 = hdulist[1].data
            for ifield in range(0, self.nfields):
                if ifield % 100 == 0:
                    print('field No. = %d/%d' % (ifield, self.nfields))
                for iBand in range(0, sdss.nBand):
                    psf = sdsspsf(hdu1, ifield, iBand)
                    psf.fit2vonK_curve_fit(vonK1arcsec)

                    # fwhmeff of the vonK1arcsec = 1.222 arcsec
                    fid.write('%d \t %d \t %d \t %5.3f \t'%(
                        ifield, camcol, iBand, psf.scaleR*1.222))
                    for i in range(len(otherParams)):
                        fid.write('%s \t'%otherParams[i][1]%(
                            hdu1[ifield][otherParams[i][0].strip()][iBand]))
                    fid.write('\n')
        fid.close()

