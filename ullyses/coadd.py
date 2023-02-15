import os
import glob

import numpy as np
from astropy.io import fits

# coadd data
#

STIS_NON_CCD_DETECTORS = ['FUV-MAMA', 'NUV-MAMA']

class SegmentList:

    def __init__(self, instrument, grating, path='.'):
        self.get_datasets = False
        if instrument is not None and grating is not None:
            self.get_datasets = True
        if grating is not None: self.grating = grating.upper()
        self.min_wavelength = None
        self.max_wavelength = None
        self.output_wavelength = None
        self.output_sumflux = None
        self.output_sumweight = None
        self.output_flux = None
        self.output_errors = None
        if instrument is not None: self.instrument = instrument.upper()
        self.detector = ''
        self.disperser = ''
        self.cenwave = ''
        self.aperture = ''
        self.s_region = ''
        self.obsmode = ''
        self.targnames = []
        self.targ_ra = ''
        self.targ_dec = ''
        self.target = ''
        self.prog_id = ''
        self.datasets = []
        self.level0 = False

        if self.get_datasets:
            vofiles = glob.glob(os.path.join(path, '*_vo.fits'))
            x1dfiles = glob.glob(os.path.join(path, '*_x1d.fits')) + glob.glob(os.path.join(path, '*_sx1.fits'))

            alldata = []
            allhdr0 = []
            allhdr1 = []

            for file in x1dfiles:
                with fits.open(file) as f1:
                    prihdr = f1[0].header
                    if prihdr['OPT_ELEM'].upper() == self.grating and prihdr['INSTRUME'].upper() == self.instrument:
                        hdr1 = f1[1].header
                        data = f1[1].data
                        alldata.append(data)
                        allhdr0.append(prihdr)
                        allhdr1.append(hdr1)
                    else:
                        continue

                if len(data) > 0:
                    print('{} added to file list for instrument/grating {}/{}'.format(file, instrument, grating))
                    self.instrument = prihdr['INSTRUME'].upper()
                    self.datasets.append(file)
                    target = prihdr['TARGNAME'].strip()
                    if target not in self.targnames:
                        self.targnames.append(target)
                    try:
                        if prihdr['HLSP_LVL'] == 0:
                            self.level0 = True
                    except:
                        pass
                else:
                    print('{} has no data'.format(file))

            if grating == 'FUSE':
                for file in vofiles:
                    with fits.open(file) as f1:
                        prihdr = f1[0].header
                        hdr1 = f1[1].header
                        data = f1[1].data
                        alldata.append(data)
                        allhdr0.append(prihdr)
                        allhdr1.append(hdr1)

                    if len(data) > 0:
                        print('{} added to file list for FUSE'.format(file))
                        self.instrument = 'FUSE'
                        aperture = prihdr["APERTURE"]
                        self.aperture = aperture
                        self.datasets.append(file)
                        target = prihdr['TARGNAME'].strip()
                        if target not in self.targnames:
                            self.targnames.append(target)
                    else:
                        print('{} has no data'.format(file))

            self.members = []
            self.primary_headers = []
            self.first_headers = []

            if len(alldata) > 0:
                for i in range(len(alldata)):
                    data = alldata[i]
                    hdr0 = allhdr0[i]
                    hdr1 = allhdr1[i]
                    if len(data) > 0:
                        self.primary_headers.append(hdr0)
                        self.first_headers.append(hdr1)
                        if self.instrument == 'FUSE':
                            sdqflags = 3
                        else:
                            sdqflags = hdr1['SDQFLAGS']
                            if self.instrument == "STIS" and (sdqflags&16) == 16 and \
                                    hdr0['DETECTOR'] in STIS_NON_CCD_DETECTORS:
                                sdqflags -= 16
                        if self.instrument == 'FUSE':
                            exptime = hdr1['EXPOSURE']
                        else:
                            exptime = hdr1['EXPTIME']
                        for row in data:
                            if self.instrument == 'COS':
                                cenwave = hdr0['CENWAVE']
                                fppos = hdr0['FPPOS']
                                if hasattr(self, "bad_segments"):
                                    if cenwave in self.bad_segments and self.bad_segments[cenwave] == row['SEGMENT']:
                                        continue

                            segment = Segment()
                            segment.data = row
                            segment.sdqflags = sdqflags
                            segment.exptime = exptime
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

        max_delta_wavelength = 0.0

        for segment in self.members:
            wavediffs = segment.data['wavelength'][1:] - segment.data['wavelength'][:-1]
            max_delta_wavelength = max(max_delta_wavelength, wavediffs.max())

        self.delta_wavelength = max_delta_wavelength

        wavegrid = np.arange(self.min_wavelength, self.max_wavelength, self.delta_wavelength)

        self.output_wavelength = wavegrid
        self.nelements = len(wavegrid)
        self.output_sumgcounts = np.zeros(self.nelements)
        self.output_sumflux = np.zeros(self.nelements)
        self.output_sumweight = np.zeros(self.nelements)
        self.output_varsum = np.zeros(self.nelements)
        self.output_flux = np.zeros(self.nelements)
        self.output_errors = np.zeros(self.nelements)
        self.signal_to_noise = np.zeros(self.nelements)
        self.sumnetcounts = np.zeros(self.nelements)
        self.output_exptime = np.zeros(self.nelements)

        return wavegrid

    def wavelength_to_index(self, wavelength):
        index = (wavelength - self.min_wavelength) / self.delta_wavelength
        indices = [int(round(x)) for x in index]
        return indices

    def index_to_wavelength(self, index):
        wavelength = index * self.delta_wavelength + self.min_wavelength
        return wavelength

    def get_flux_weight(self, segment):
        pass

    def coadd(self):
        for segment in self.members:
            goodpixels = np.where((segment.data['dq'] & segment.sdqflags) == 0)
            wavelength = segment.data['wavelength'][goodpixels]
            indices = self.wavelength_to_index(wavelength)
            npts = len(indices)
            flux_weight = self.get_flux_weight(segment)
            gcounts0 = np.abs(segment.data['gross'] * segment.exptime)
            gcounts = gcounts0[goodpixels]
            weight = flux_weight[goodpixels]
            net_counts = segment.data['net'][goodpixels] * segment.exptime
            if self.instrument == 'COS':
                variances = segment.data['variance_counts'][goodpixels] + segment.data['variance_bkg'][goodpixels] + segment.data['variance_flat'][goodpixels]

            flux = segment.data['flux'][goodpixels]
            for i in range(npts):
                self.sumnetcounts[indices[i]] = self.sumnetcounts[indices[i]] + net_counts[i]
                self.output_sumweight[indices[i]] = self.output_sumweight[indices[i]] + weight[i]
                self.output_sumgcounts[indices[i]] = self.output_sumgcounts[indices[i]] + gcounts[i]
                self.output_sumflux[indices[i]] = self.output_sumflux[indices[i]] + flux[i] * weight[i]
                self.output_exptime[indices[i]] = self.output_exptime[indices[i]] + segment.exptime
                if self.instrument == 'COS':
                    self.output_varsum[indices[i]] = self.output_varsum[indices[i]] + variances[i]
        good_dq = np.where(self.output_exptime > 0.)
        self.first_good_wavelength = self.output_wavelength[good_dq][0]
        self.last_good_wavelength = self.output_wavelength[good_dq][-1]
        nonzeros = np.where(self.output_sumweight != 0)
        if self.instrument == 'COS':
            # Using the variances (which only COS has) gives spikes in the error when the flux goes negative.
            self.output_varsum[nonzeros] = np.where(self.output_varsum[nonzeros] < 0.5, 0.5, self.output_varsum[nonzeros])
            self.output_errors[nonzeros] = np.sqrt(self.output_varsum[nonzeros])
        else:
            # For the moment calculate errors from the gross counts
            self.output_errors[nonzeros] = np.sqrt(self.output_sumgcounts[nonzeros])
        self.output_flux[nonzeros] = self.output_sumflux[nonzeros] / self.output_sumweight[nonzeros]
        # Calculate the conversion from counts to flux units for errors
        conversion = self.output_flux[nonzeros] / self.sumnetcounts[nonzeros]
        # Clean out NaNs from where flux and net are zero
        good = np.where(~np.isnan(conversion))
        bad = np.where(np.isnan(conversion))
        lenwl = len(wavelength)
        interpolated_values = np.interp(self.output_wavelength[bad],
                                        self.output_wavelength[good],
                                        conversion[good])
        conversion[bad] = interpolated_values
        # Use conversion to calculate error in flux units (erg/cm^2/s/Ang)
        self.output_errors[nonzeros] = self.output_errors[nonzeros] * conversion
        self.signal_to_noise[nonzeros] = self.output_flux[nonzeros] / self.output_errors[nonzeros]
        return

    def get_targname(self):
        # These are just preliminary target names, in case we can't find a match
        if len(self.targnames) == 1:
            targname = self.targnames[0]
        else:
            targname = self.primary_headers[0]["targname"]
        return targname

    def get_coords(self):
        ras = list(set([h["ra_targ"] for h in self.primary_headers]))
        decs = list(set([h["dec_targ"] for h in self.primary_headers]))
        avg_ra = np.average(ras)
        avg_dec = np.average(decs)
        return avg_ra, avg_dec

# Weight functions for STIS
weight_function = {
    'unity':      lambda x, y, z: np.ones(len(x)),
    'gross':      lambda x, y, z: x,
    'exptime':    lambda x, y, z: x * y,
    'throughput': lambda x, y, z: y * z
}

class STISSegmentList(SegmentList):

    def __init__(self, grating, path, weighting_method='throughput'):
        instrument = 'STIS'
        if grating is None: instrument = None
        SegmentList.__init__(self, instrument, grating, path)

        self.weighting_method = weighting_method

    def get_flux_weight(self, segment):
       exptime = segment.exptime
       gross = segment.data['gross']
       net = segment.data['net']
       flux = segment.data['flux']

       weighted_gross = weight_function[self.weighting_method](gross, exptime, net/flux)

       return np.abs(weighted_gross)


class COSSegmentList(SegmentList):

    def __init__(self, grating, path):
        instrument = 'COS'
        if grating is None: instrument = None
        SegmentList.__init__(self, instrument, grating, path)

    def get_flux_weight(self, segment):
        thru_nans = segment.data['net'] / segment.data['flux']
        if set(np.isnan(thru_nans)) == {False}:
            weight = thru_nans * segment.exptime
        else:
            xpix = np.arange(len(thru_nans))
            good = np.where( np.isnan(thru_nans) == False)
            good_xpix = xpix[good]
            good_thru_nans = thru_nans[good]
            thru = np.interp(xpix, good_xpix, good_thru_nans)
            weight = thru * segment.exptime

        return np.abs(weight)

class FUSESegmentList(SegmentList):

    def __init__(self, grating, path):
        SegmentList.__init__(self, 'FUSE', 'FUSE', path)

    def get_gross_counts(self, segment):
        print("FUSE doesn't use gross counts")
        return

    def create_output_wavelength_grid(self):
        #
        # For FUSE data, we just need to return the wavelength grid associated
        # with the data, and make sure we appropriately populate the same attributes
        # as the base class
        segment = self.members[0]
        self.min_wavelength = segment.data['wave'].min()
        self.max_wavelength = segment.data['wave'].max()

        self.delta_wavelength = None

        self.output_wavelength = segment.data['wave']
        self.nelements = len(self.output_wavelength)
        self.output_sumflux = np.zeros(self.nelements)
        self.output_sumweight = np.zeros(self.nelements)
        self.output_flux = np.zeros(self.nelements)
        self.output_errors = np.zeros(self.nelements)
        self.signal_to_noise = np.zeros(self.nelements)
        self.output_exptime = np.zeros(self.nelements)

        return self.output_wavelength

    def coadd(self):

        segment = self.members[0]
        nelements = len(self.output_wavelength)
        try:
            dq = segment.data['dq']
        except KeyError:
            dq = np.ones(nelements).astype(np.int32)
        goodpixels = np.where((dq & segment.sdqflags) == 0)
        self.output_exptime[goodpixels] = segment.exptime
        self.output_flux[goodpixels] = segment.data['flux'][goodpixels]
        self.output_errors[goodpixels] = segment.data['sigma'][goodpixels]
        nonzeros = np.where(self.output_errors != 0.0)
        self.signal_to_noise[nonzeros] = self.output_flux[nonzeros] / self.output_errors[nonzeros]
        good_dq = np.where(self.output_exptime > 0.)
        self.first_good_wavelength = self.output_wavelength[good_dq][0]
        self.last_good_wavelength = self.output_wavelength[good_dq][-1]
        return

class CCDSegmentList(SegmentList):

    def __init__(self, grating, path, weighting_method='throughput'):
        instrument = 'STIS'
        if grating is None: instrument = None
        SegmentList.__init__(self, instrument, grating, path)

        self.weighting_method = weighting_method

    def get_flux_weight(self, segment):
       exptime = segment.exptime
       gross = segment.data['gross']
       net = segment.data['net']
       flux = segment.data['flux']

       weighted_gross = weight_function[self.weighting_method](gross, exptime, net/flux)

       return np.abs(weighted_gross)

    def create_output_wavelength_grid(self):
        if len(self.members) == 1:
            #
            # For STIS CCD data when there is only 1 segment, we just need
            # to return the wavelength grid associated
            # with the data, and make sure we appropriately populate the same attributes
            # as the base class
            segment = self.members[0]
            self.min_wavelength = segment.data['wavelength'].min()
            self.max_wavelength = segment.data['wavelength'].max()

            self.delta_wavelength = None

            self.output_wavelength = segment.data['wavelength']
        else:
            min_wavelength = 20000.0
            max_wavelength = 0.0
            for segment in self.members:
                minwave = segment.data['wavelength'].min()
                maxwave = segment.data['wavelength'].max()
                if minwave < min_wavelength: min_wavelength = minwave
                if maxwave > max_wavelength: max_wavelength = maxwave
            self.min_wavelength = int(min_wavelength)
            self.max_wavelength = int(max_wavelength) + 1

            max_delta_wavelength = 0.0

            for segment in self.members:
                wavediffs = segment.data['wavelength'][1:] - segment.data['wavelength'][:-1]
                max_delta_wavelength = max(max_delta_wavelength, wavediffs.max())

            self.delta_wavelength = max_delta_wavelength

            wavegrid = np.arange(self.min_wavelength, self.max_wavelength, self.delta_wavelength)

            self.output_wavelength = wavegrid

        self.nelements = len(self.output_wavelength)
        self.output_sumflux = np.zeros(self.nelements)
        self.output_sumerrors = np.zeros(self.nelements)
        self.output_sumweight = np.zeros(self.nelements)
        self.output_flux = np.zeros(self.nelements)
        self.output_errors = np.zeros(self.nelements)
        self.signal_to_noise = np.zeros(self.nelements)
        self.output_exptime = np.zeros(self.nelements)

        return self.output_wavelength

    def coadd(self):
        if len(self.members) == 1:
            segment = self.members[0]
            nelements = len(self.output_wavelength)
            try:
                dq = segment.data['dq']
            except KeyError:
                dq = np.ones(nelements).astype(np.int32)
            goodpixels = np.where((dq & segment.sdqflags) == 0)
            self.output_exptime[goodpixels] = segment.exptime
            self.output_flux[goodpixels] = segment.data['flux'][goodpixels]
            self.output_errors[goodpixels] = segment.data['error'][goodpixels]
            nonzeros = np.where(self.output_errors != 0.0)
            self.signal_to_noise[nonzeros] = self.output_flux[nonzeros] / self.output_errors[nonzeros]
            good_dq = np.where(self.output_exptime > 0.)
            self.first_good_wavelength = self.output_wavelength[good_dq][0]
            self.last_good_wavelength = self.output_wavelength[good_dq][-1]
        else:
            nelements = len(self.output_wavelength)
            self.output_sumsqerrors = np.zeros(nelements)
            for segment in self.members:
                goodpixels = np.where((segment.data['dq'] & segment.sdqflags) == 0)
                wavelength = segment.data['wavelength'][goodpixels]
                indices = self.wavelength_to_index(wavelength)
                npts = len(indices)
                weight = self.get_flux_weight(segment)[goodpixels]
                flux = segment.data['flux'][goodpixels]
                err = segment.data['error'][goodpixels]
                for i in range(npts):
                    self.output_sumweight[indices[i]] = self.output_sumweight[indices[i]] + weight[i]
                    self.output_sumflux[indices[i]] = self.output_sumflux[indices[i]] + flux[i] * weight[i]
                    self.output_sumsqerrors[indices[i]] = self.output_sumsqerrors[indices[i]] + (err[i] * weight[i])**2
                    self.output_exptime[indices[i]] = self.output_exptime[indices[i]] + segment.exptime
            good_dq = np.where(self.output_exptime > 0.)
            self.first_good_wavelength = self.output_wavelength[good_dq][0]
            self.last_good_wavelength = self.output_wavelength[good_dq][-1]
            nonzeros = np.where(self.output_sumweight != 0)
            self.output_flux[nonzeros] = self.output_sumflux[nonzeros] / self.output_sumweight[nonzeros]
            self.output_errors[nonzeros] = np.sqrt(self.output_sumsqerrors[nonzeros]) / self.output_sumweight[nonzeros]
            self.signal_to_noise[nonzeros] = self.output_flux[nonzeros] / self.output_errors[nonzeros]
#            self.output_errors[nonzeros] = np.abs(self.output_flux[nonzeros] / self.signal_to_noise[nonzeros])
        return

    def get_inverse_variance(self, segment):
        error = segment.data['error']
        inverse_variance = 1.0 / (error*error)
        return inverse_variance

class Segment:

    def __init__(self):
        self.data = None
        self.sdqflags = None
        self.exptime = None


def abut(product_short, product_long):
    """Abut the spectra in 2 products.  Assumes the first argument is the shorter
    wavelength.  If either product is None, just return the product that isn't None.
    """
    if product_short is not None and product_long is not None:
        transition_wavelength = find_transition_wavelength(product_short, product_long)
        if transition_wavelength == "bad":
            return None
        # Spectra are overlapped
        if transition_wavelength is not None:
            short_indices = np.where(product_short.output_wavelength < transition_wavelength)
            transition_index_short = short_indices[0][-1]
            long_indices = np.where(product_long.output_wavelength > transition_wavelength)
            transition_index_long = long_indices[0][0]
        else:
            # No overlap
            goodshort = np.where(product_short.output_exptime > 0.)
            goodlong = np.where(product_long.output_exptime > 0.)
            last_good_short = product_short.output_wavelength[goodshort][-1]
            first_good_long = product_long.output_wavelength[goodlong][0]
            short_indices = np.where(product_short.output_wavelength < last_good_short)
            transition_index_short = short_indices[0][-1]
            long_indices = np.where(product_long.output_wavelength > first_good_long)
            transition_index_long = long_indices[0][0]
        output_grating = product_short.grating + '-' + product_long.grating
        product_abutted = product_short.__class__('', output_grating)
        nout = len(product_short.output_wavelength[:transition_index_short])
        nout = nout + len(product_long.output_wavelength[transition_index_long:])
        product_abutted.nelements = nout
        product_abutted.output_wavelength = np.zeros(nout)
        product_abutted.output_flux = np.zeros(nout)
        product_abutted.output_errors = np.zeros(nout)
        product_abutted.signal_to_noise = np.zeros(nout)
        product_abutted.output_exptime = np.zeros(nout)
        product_abutted.output_wavelength[:transition_index_short] = product_short.output_wavelength[:transition_index_short]
        product_abutted.output_wavelength[transition_index_short:] = product_long.output_wavelength[transition_index_long:]
        product_abutted.output_flux[:transition_index_short] = product_short.output_flux[:transition_index_short]
        product_abutted.output_flux[transition_index_short:] = product_long.output_flux[transition_index_long:]
        product_abutted.output_errors[:transition_index_short] = np.abs(product_short.output_errors[:transition_index_short])
        product_abutted.output_errors[transition_index_short:] = np.abs(product_long.output_errors[transition_index_long:])
        product_abutted.signal_to_noise[:transition_index_short] = product_short.signal_to_noise[:transition_index_short]
        product_abutted.signal_to_noise[transition_index_short:] = product_long.signal_to_noise[transition_index_long:]
        product_abutted.output_exptime[:transition_index_short] = product_short.output_exptime[:transition_index_short]
        product_abutted.output_exptime[transition_index_short:] = product_long.output_exptime[transition_index_long:]
        product_abutted.primary_headers = product_short.primary_headers + product_long.primary_headers
        product_abutted.first_headers = product_short.first_headers + product_long.first_headers
        product_abutted.grating = output_grating
        product_short.target = product_short.get_targname()
        product_long.target = product_long.get_targname()
        if product_short.instrument in product_long.instrument:
            product_abutted.instrument = product_long.instrument
        if product_long.instrument in product_short.instrument:
            product_abutted.instrument = product_short.instrument
        else:
            product_abutted.instrument = product_short.instrument + '-' + product_long.instrument
        target_matched = False
        if product_short.target == product_long.target:
            product_abutted.target = product_short.target
            target_matched = True
            product_abutted.targnames = [product_short.target, product_long.target]
        else:
            for target_name in product_short.targnames:
                if target_name in product_long.targnames:
                    product_abutted.target = product_abutted.get_targname()
                    target_matched = True
                    product_abutted.targnames = [target_name]
        if not target_matched:
            product_abutted = None
            print(f'Trying to abut spectra from 2 different targets:')
            print(f'{product_short.target} and {product_long.target}')
    else:
        if product_short is not None:
            product_abutted = product_short
        elif product_long is not None:
            product_abutted = product_long
        else:
            product_abutted = None
    product_abutted.targ_ra, product_abutted.targ_dec = product_abutted.get_coords()
    return product_abutted


def find_transition_wavelength(product_short, product_long):
    """Find the wavelength below which we use product_short and above
    which we use product_long.
    Initial implementation is to return the wavelength midway between the last good short wavelength and the first
    good long wavelength.  If there's no overlap, return None
    If the long product is totally encompassed in the short product's
    wavelength range, return 'bad' so no abutting is attempted.
    """

    goodshort = np.where(product_short.output_exptime > 0.)
    goodlong = np.where(product_long.output_exptime > 0.)
    first_good_short = product_short.output_wavelength[goodshort][0]
    last_good_short = product_short.output_wavelength[goodshort][-1]
    first_good_long = product_long.output_wavelength[goodlong][0]
    last_good_long = product_long.output_wavelength[goodlong][-1]
    if last_good_long < last_good_short: # long product is entirely in short spectrum
        return "bad"
    if first_good_short > first_good_long and last_good_short < last_good_long: # short is entirely in long
        return "bad"
    if last_good_short > first_good_long:
        return 0.5*(last_good_short + first_good_long)
    else:
        return None
