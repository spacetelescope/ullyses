import os
import glob
import datetime
from datetime import datetime as dt
import re

import pandas as pd
import numpy as np
import astropy
from astropy.io import fits
from astropy.time import Time

import ullyses_utils
# coadd data
#

CAL_VER = 1.1
RED = "\033[1;31m"
RESET = "\033[0;0m"


STIS_NON_CCD_DETECTORS = ['FUV-MAMA', 'NUV-MAMA']
BAD_SEGMENTS = ['NUVC']

class SegmentList:

    def __init__(self, instrument, grating, path='.'):
        self.grating = grating.upper()
        self.min_wavelength = None
        self.max_wavelength = None
        self.output_wavelength = None
        self.output_sumflux = None
        self.output_sumweight = None
        self.output_flux = None
        self.output_errors = None
        self.instrument = instrument.upper()
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
                        if self.instrument == 'COS' and row['SEGMENT'] in BAD_SEGMENTS:
                            continue

                        segment = Segment()
                        segment.data = row
                        segment.sdqflags = sdqflags
                        segment.exptime = exptime
                        if self.instrument == 'COS':
                            cenwave = hdr0['CENWAVE']
                            fppos = hdr0['FPPOS']
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

    def write(self, filename, overwrite=False, level="", version=""):

        # If the target is a ULLYSES target, use the official
        # target name and coords
        self.target = self.ull_targname()
        self.targ_ra, self.targ_dec = self.ull_coords()

        # Table 1 - HLSP data

        # set up the header
        hdr1 = fits.Header()
        hdr1['EXTNAME'] = ('SCIENCE', 'Spectrum science arrays')
        hdr1['TIMESYS'] = ('UTC', 'Time system in use')
        hdr1['TIMEUNIT'] = ('s', 'Time unit for durations')
        hdr1['TREFPOS'] = ('GEOCENTER', 'Time reference position')

        mjd_beg = self.combine_keys("expstart", "min")
        mjd_end = self.combine_keys("expend", "max")
        dt_beg = Time(mjd_beg, format="mjd").datetime
        dt_end = Time(mjd_end, format="mjd").datetime
        hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S"), 'Date-time of first observation start')
        hdr1.add_blank('', after='TREFPOS')
        hdr1.add_blank('              / FITS TIME COORDINATE KEYWORDS', before='DATE-BEG')

        hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
        hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
        hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
        hdr1['XPOSURE'] = (self.combine_keys("exptime", "sum"), '[s] Sum of exposure durations')

        # set up the table columns
        nelements = len(self.output_wavelength)
        rpt = str(nelements)

        # Table with co-added spectrum
        cw = fits.Column(name='WAVELENGTH', format=rpt+'E', unit="Angstrom")
        cf = fits.Column(name='FLUX', format=rpt+'E', unit="erg /s /cm**2 /Angstrom")
        ce = fits.Column(name='ERROR', format=rpt+'E', unit="erg /s /cm**2 /Angstrom")
        cs = fits.Column(name='SNR', format=rpt+'E')
        ct = fits.Column(name='EFF_EXPTIME', format=rpt+'E', unit="s")
        cd = fits.ColDefs([cw, cf, ce, cs, ct])
        table1 = fits.BinTableHDU.from_columns(cd, nrows=1, header=hdr1)

        # populate the table
        table1.data['WAVELENGTH'] = self.output_wavelength.copy()
        table1.data['FLUX'] = self.output_flux.copy()
        table1.data['ERROR'] = self.output_errors.copy()
        table1.data['SNR'] = self.signal_to_noise.copy()
        table1.data['EFF_EXPTIME'] = self.output_exptime.copy()
        # HLSP primary header
        hdr0 = fits.Header()
        hdr0['EXTEND'] = ('T', 'FITS file may contain extensions')
        hdr0['NEXTEND'] = 3
        hdr0['FITS_VER'] = 'Definition of the Flexible Image Transport System (FITS) v4.0 https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf'
        hdr0['FITS_SW'] = ('astropy.io.fits v' + astropy.__version__, 'FITS file creation software')
        hdr0['ORIGIN'] = ('Space Telescope Science Institute', 'FITS file originator')
        hdr0['DATE'] = (str(datetime.date.today()), 'Date this file was written')
        hdr0['FILENAME'] = (os.path.basename(filename), 'Name of this file')
        hdr0['TELESCOP'] = (self.combine_keys("telescop", "multi"), 'Telescope used to acquire data')
        hdr0['INSTRUME'] = (self.combine_keys("instrume", "multi"), 'Instrument used to acquire data')
        hdr0.add_blank('', after='TELESCOP')
        hdr0.add_blank('              / SCIENCE INSTRUMENT CONFIGURATION', before='INSTRUME')
        hdr0['DETECTOR'] = (self.combine_keys("detector", "multi"), 'Detector or channel used to acquire data')
        hdr0['DISPERSR'] = (self.combine_keys("opt_elem", "multi"), 'Identifier of disperser')
        hdr0['CENWAVE'] = (self.combine_keys("cenwave", "multi"), 'Central wavelength setting for disperser')
        hdr0['APERTURE'] = (self.combine_keys("aperture", "multi"), 'Identifier of entrance aperture')
        hdr0['S_REGION'] = (self.obs_footprint(), 'Region footprint')
        hdr0['OBSMODE'] = (self.combine_keys("obsmode", "multi"), 'Instrument operating mode (ACCUM | TIME-TAG)')
        hdr0['TARGNAME'] = self.target
        hdr0.add_blank(after='OBSMODE')
        hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

        hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
        hdr0['TARG_RA'] =  (self.targ_ra,  '[deg] Target right ascension')
        hdr0['TARG_DEC'] =  (self.targ_dec,  '[deg] Target declination')
        hdr0['PROPOSID'] = (self.combine_keys("proposid", "multi"), 'Program identifier')
        hdr0.add_blank(after='TARG_DEC')
        hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
        hdr0['CAL_VER'] = (f'ULLYSES Cal {CAL_VER}', 'HLSP processing software version')
        hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
        hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                        'Name ID of this HLSP collection')
        hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead')
        hdr0['HLSP_VER'] = (version,'HLSP data release version identifier')
        hdr0['HLSP_LVL'] = (level, 'ULLYSES HLSP Level')
        hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
        hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
        hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')

        hdr0['CENTRWV'] = (self.combine_keys("centrwv", "average"), 'Central wavelength of the data')
        hdr0.add_blank(after='REFERENC')
        hdr0.add_blank('           / ARCHIVE SEARCH KEYWORDS', before='CENTRWV')
        hdr0['MINWAVE'] = (self.combine_keys("minwave", "min"), 'Minimum wavelength in spectrum')
        hdr0['MAXWAVE'] = (self.combine_keys("maxwave", "max"), 'Maximum wavelength in spectrum')

        primary = fits.PrimaryHDU(header=hdr0)

        # Table 2 - individual product info

        # first set up header
        hdr2 = fits.Header()
        hdr2['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
        # set up the table columns
        cfn = fits.Column(name='FILENAME', array=self.combine_keys("filename", "arr"), format='A64')
        cpid = fits.Column(name='PROPOSID', array=self.combine_keys("proposid", "arr"), format='A32')
        ctel = fits.Column(name='TELESCOPE', array=self.combine_keys("telescop", "arr"), format='A32')
        cins = fits.Column(name='INSTRUMENT', array=self.combine_keys("instrume", "arr"), format='A32')
        cdet = fits.Column(name='DETECTOR', array=self.combine_keys("detector", "arr"), format='A32')
        cdis = fits.Column(name='DISPERSER', array=self.combine_keys("opt_elem", "arr"), format='A32')
        ccen = fits.Column(name='CENWAVE', array=self.combine_keys("cenwave", "arr"), format='A32')
        cap = fits.Column(name='APERTURE', array=self.combine_keys("aperture", "arr"), format='A32')
        csr = fits.Column(name='SPECRES', array=self.combine_keys("specres", "arr"), format='F8.1')
        ccv = fits.Column(name='CAL_VER', array=self.combine_keys("cal_ver", "arr"), format='A32')
        mjd_begs = self.combine_keys("expstart", "arr")
        mjd_ends = self.combine_keys("expend", "arr")
        mjd_mids = (mjd_ends + mjd_begs) / 2.
        cdb = fits.Column(name='MJD_BEG', array=mjd_begs, format='F15.9', unit='d')
        cdm = fits.Column(name='MJD_MID', array=mjd_mids, format='F15.9', unit='d')
        cde = fits.Column(name='MJD_END', array=mjd_ends, format='F15.9', unit='d')
        cexp = fits.Column(name='XPOSURE', array=self.combine_keys("exptime", "arr"), format='F15.9', unit='s')
        cmin = fits.Column(name='MINWAVE', array=self.combine_keys("minwave", "arr"), format='F9.4', unit='Angstrom')
        cmax = fits.Column(name='MAXWAVE', array=self.combine_keys("maxwave", "arr"), format='F9.4', unit='Angstrom')

        cd2 = fits.ColDefs([cfn, cpid, ctel, cins, cdet, cdis, ccen, cap, csr, ccv, cdb, cdm, cde, cexp, cmin ,cmax])

        table2 = fits.BinTableHDU.from_columns(cd2, header=hdr2)

        # the HDUList:
        # 0 - empty data - 0th ext header
        # 1 - HLSP data - 1st ext header
        # 2 - individual product info - 2nd ext header

        hdul = fits.HDUList([primary, table1, table2])

        hdul.writeto(filename, overwrite=overwrite)

        # from ullyses_jira.parse_csv import parse_name_csv
        # name_mapping = {}
        # for ttype in ['lmc', 'smc', 'tts']:
        #     names_dict = parse_name_csv(ttype)
        #     name_mapping = {**name_mapping, **names_dict}


    def obs_footprint(self):
        # Not using WCS at the moment
        # This is a placeholder, need to figure out polygon
#        apertures = list(set([h["aperture"] for h in self.primary_headers]))
#        ras = list(set([h["ra_targ"] for h in self.primary_headers]))
#        ra_diff = max(np.abs(ras)) - min(np.abs(ras))
#        decs = list(set([h["dec_targ"] for h in self.primary_headers]))
#        dec_diff = max(np.abs(decs)) - min(np.abs(decs))
#        center_ra = np.average(ras)
#        center_dec = np.average(decs)
#        extent_ra = (2.5 / 2 / 3600) + ra_diff
#        extent_dec = (2.5 / 2 / 3600) + dec_diff
#        radius = max([extent_ra, extent_dec])
        self.targ_ra, self.targ_dec = self.ull_coords()
        radius = (2.5 / 2 / 3600)
        center_ra = self.targ_ra
        center_dec = self.targ_dec

        s_region = f"CIRCLE {center_ra} {center_dec} {radius}"
        return s_region

    def ull_targname(self):
        aliases_file = ullyses_utils.__path__[0] + '/data/target_metadata/pd_all_aliases.json'
        aliases = pd.read_json(aliases_file, orient="split")
        # These are just preliminary target names, in case we can't find a match
        if len(self.targnames) == 1:
            ull_targname = self.targnames[0]
        else:
            ull_targname = self.primary_headers[0]["targname"]
        targ_matched = False
        for targ in self.targnames:
            # The alias file is all uppercase
            targ_upper = targ.upper()
            mask = aliases.apply(lambda row: row.astype(str).str.fullmatch(re.escape(targ_upper)).any(), axis=1)
            if set(mask) != {False}:
                targ_matched = True
                ull_targname = aliases[mask]["ULL_MAST_name"].values[0]
                break
        if targ_matched is False:
            print(f"{RED}WARNING: Could not match target name {ull_targname} to ULLYSES alias list{RESET}")
        return ull_targname

    def ull_coords(self):
        ras = list(set([h["ra_targ"] for h in self.primary_headers]))
        decs = list(set([h["dec_targ"] for h in self.primary_headers]))
        avg_ra = np.average(ras)
        avg_dec = np.average(decs)
        if self.target == "":
            return avg_ra, avg_dec
        targetinfo_file = ullyses_utils.__path__[0] + '/data/target_metadata/pd_targetinfo.json'
        master_list = pd.read_json(targetinfo_file, orient="split")
        master_list = master_list.apply(lambda x: x.astype(str).str.upper())
        coords = master_list.loc[master_list["mast_targname"] == self.target][["ra", "dec"]].values
        if len(coords) != 0:
            return coords[0][0], coords[0][1]
        else:
            return avg_ra, avg_dec


    def combine_keys(self, key, method):
        keymap= {"HST": {"expstart": ("expstart", 1),
                         "expend": ("expend", 1),
                         "exptime": ("exptime", 1),
                         "telescop": ("telescop", 0),
                         "instrume": ("instrume", 0),
                         "detector": ("detector", 0),
                         "opt_elem": ("opt_elem", 0),
                         "cenwave": ("cenwave", 0),
                         "aperture": ("aperture", 0),
                         "obsmode": ("obsmode", 0),
                         "proposid": ("proposid", 0),
                         "centrwv": ("centrwv", 0),
                         "minwave": ("minwave", 0),
                         "maxwave": ("maxwave", 0),
                         "filename": ("filename", 0),
                         "specres": ("specres", 0),
                         "cal_ver": ("cal_ver", 0)},
                "FUSE": {"expstart": ("obsstart", 0),
                         "expend": ("obsend", 0),
                         "exptime": ("obstime", 0),
                         "telescop": ("telescop", 0),
                         "instrume": ("instrume", 0),
                         "detector": ("detector", 0),
                         "opt_elem": ("detector", 0),
                         "cenwave": ("centrwv", 0),
                         "aperture": ("aperture", 0),
                         "obsmode": ("instmode", 0),
                         "proposid": ("prgrm_id", 0),
                         "centrwv": ("centrwv", 0),
                         "minwave": ("wavemin", 0),
                         "maxwave": ("wavemax", 0),
                         "filename": ("filename", 0),
                         "specres": ("spec_rp", 1),
                         "cal_ver": ("cf_vers", 0)}}

        vals = []
        for i in range(len(self.primary_headers)):
            tel = self.primary_headers[i]["telescop"]
            actual_key = keymap[tel][key][0]
            hdrno = keymap[tel][key][1]
            if hdrno == 0:
                val = self.primary_headers[i][actual_key]
            else:
                val = self.first_headers[i][actual_key]
            if tel == "FUSE" and key == "filename":
                val = val.replace(".fit", "_vo.fits")
            vals.append(val)

        # Allowable methods are min, max, average, sum, multi, arr
        if method == "multi":
            keys_set = list(set(vals))
            if len(keys_set) > 1:
                return "MULTI"
            else:
                return keys_set[0]
        elif method == "min":
            return min(vals)
        elif method == "max":
            return max(vals)
        elif method == "average":
            return np.average(vals)
        elif method == "sum":
            return np.sum(vals)
        elif method == "arr":
            return np.array(vals)


# Weight functions for STIS
weight_function = {
    'unity':      lambda x, y, z: np.ones(len(x)),
    'gross':      lambda x, y, z: x,
    'exptime':    lambda x, y, z: x * y,
    'throughput': lambda x, y, z: y * z
}

class STISSegmentList(SegmentList):

    def __init__(self, grating, path, weighting_method='throughput'):
        SegmentList.__init__(self, 'STIS', grating, path)

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
        SegmentList.__init__(self, 'COS', grating, path)

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
        SegmentList.__init__(self, 'STIS', grating, path)

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
        product_abutted = SegmentList('', output_grating)
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
        product_short.target = product_short.ull_targname()
        product_long.target = product_long.ull_targname()
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
                    product_abutted.target = product_abutted.ull_targname()
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
    product_abutted.targ_ra, product_abutted.targ_dec = product_abutted.ull_coords()
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
