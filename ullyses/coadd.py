import os
import glob

import numpy as np

from astropy.io import fits

import parameters

#
# coadd data
#


class COSSegmentList:

    def __init__(self, grating, path='.'):
        self.grating = grating
        self.min_wavelength = None
        self.max_wavelength = None
        self.output_wavelength = None
        self.output_sumflux = None
        self.output_sumweight = None
        self.output_flux = None
        
        x1dfiles = glob.glob(os.path.join(path, '*_x1d.fits'))

        gratinglist = []

        for file in x1dfiles:
            f1 = fits.open(file)
            prihdr = f1[0].header
            if prihdr['OPT_ELEM'] == grating:
                gratinglist.append(f1)
            else:
                f1.close()

        self.members = []
        self.primary_headers = []

        for hdulist in gratinglist:

            self.primary_headers.append(hdulist[0].header)

            data = hdulist[1].data
            sdqflags = hdulist[1].header['SDQFLAGS']
            for row in data:
                segment = Segment()
                segment.data = row
                segment.sdqflags = sdqflags
                self.members.append(segment)

    def create_output_wavelength_grid(self):
        min_wavelength = 10000.0
        max_wavelength = 0.0
        for segment in self.members:
            minwave = segment.data['wavelength'].min()
            maxwave = segment.data['wavelength'].max()
            if minwave < min_wavelength: min_wavelength = minwave
            if maxwave > max_wavelength: max_wavelength = maxwave
        self.min_wavelength = int(min_wavelength)
        self.max_wavelength = int(max_wavelength) + 1
    
        deltasum = 0.0
    
        for segment in self.members:
            wavediffs = segment.data['wavelength'][1:] - segment.data['wavelength'][:-1]
            deltasum += wavediffs.mean()
    
        self.delta_wavelength = deltasum / len(self.members)
    
        wavegrid = np.arange(self.min_wavelength, self.max_wavelength, self.delta_wavelength)
    
        self.output_wavelength = wavegrid
        self.nelements = len(wavegrid)
        self.output_sumflux = np.zeros(self.nelements)
        self.output_sumweight = np.zeros(self.nelements)
        self.output_flux = np.zeros(self.nelements)

        return wavegrid

    def wavelength_to_index(self, wavelength):
        index = (wavelength - self.min_wavelength) / self.delta_wavelength
        return index.astype(np.int)

    def index_to_wavelength(self, index):
        wavelength = index * self.delta_wavelength + self.min_wavelength
        return wavelength

    def coadd(self):
        for segment in self.members:
            goodpixels = np.where((segment.data['dq'] & segment.sdqflags) == 0)
            wavelength = segment.data['wavelength'][goodpixels]
            indices = self.wavelength_to_index(wavelength)
            weight = segment.data['gcounts'][goodpixels]
            flux = segment.data['flux'][goodpixels]
            self.output_sumweight[indices] = self.output_sumweight[indices] + weight
            self.output_sumflux[indices] = self.output_sumflux[indices] + flux * weight
        
        nonzeros = np.where(self.output_sumweight != 0)
        self.output_flux[nonzeros] = self.output_sumflux[nonzeros] / self.output_sumweight[nonzeros]
        return

    def write(self, filename, overwrite=False):

        # Table with co-added spectrum
        cw = fits.Column(name='WAVELENGTH', array=self.output_wavelength, format='E')
        cf = fits.Column(name='FLUX', array=self.output_flux, format='E')
        table = fits.BinTableHDU.from_columns([cw, cf])

        # HLSP primary header
        hdr = fits.Header()
        hdr['COMMENT'] = "Mock HLSP file."
        hdr['NEXTEND'] = len(self.primary_headers) + 1
        primary = fits.PrimaryHDU(header=hdr)

        # HLSP file is comprised of a list of HDUs, with only the first
        # one being used to store the spectrum. Remaining HDUs contain
        # the primary headers from each input spectrum. Their data sections
        # are empty. (this will likely be needed to populate the quicklook
        # tool reporting widgets)
        hdul = fits.HDUList([primary, table])

        for p_header in self.primary_headers:
            extension = fits.BinTableHDU(header=p_header)
            hdul.append(extension)

        hdul.writeto(filename, overwrite=overwrite)

class Segment:

    def __init__(self):
        self.data = None
        self.sdqflags = None

