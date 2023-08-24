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

#
# coadd data
#

CAL_VER = 0.1
RED = "\033[1;31m"
RESET = "\033[0;0m"


class SegmentList:

    def __init__(self, grating, path='.'):
        self.grating = grating
        self.min_wavelength = None
        self.max_wavelength = None
        self.output_wavelength = None
        self.output_sumflux = None
        self.output_sumweight = None
        self.output_flux = None
        self.output_errors = None
        self.instrument = None
        self.detector = ''
        self.disperser = ''
        self.cenwave = ''
        self.aperture = ''
        self.s_region = ''
        self.obsmode = ''
        self.targname = []
        self.targ_ra = ''
        self.targ_dec = ''
        self.target = ''
        self.prog_id = ''
        self.datasets = []

        x1dfiles = glob.glob(os.path.join(path, '*_x1d.fits'))

        gratinglist = []

        for file in x1dfiles:
            f1 = fits.open(file)
            prihdr = f1[0].header
            if prihdr['OPT_ELEM'] == grating:
                data = f1[1].data
                if len(data) > 0:
                    print('{} added to file list for grating {}'.format(file, grating))
                    gratinglist.append(f1)
                    self.instrument = prihdr['INSTRUME']
                    self.datasets.append(file)
                    target = prihdr['TARGNAME']
                    if target not in self.targname:
                        self.targname.append(target.strip())
                else:
                    print('{} has no data'.format(file))
            else:
                f1.close()

        self.members = []
        self.primary_headers = []
        self.first_headers = []

        if len(gratinglist) > 0:
            for hdulist in gratinglist:

                data = hdulist[1].data
                if len(data) > 0:
                    self.primary_headers.append(hdulist[0].header)
                    self.first_headers.append(hdulist[1].header)
                    sdqflags = hdulist[1].header['SDQFLAGS']
                    if self.instrument == "STIS" and (sdqflags&16) == 16:
                        sdqflags -= 16
                    exptime = hdulist[1].header['EXPTIME']
                    for row in data:
                        segment = Segment()
                        segment.data = row
                        segment.sdqflags = sdqflags
                        segment.exptime = exptime
                        if self.instrument == 'COS':
                            cenwave = hdulist[0].header['CENWAVE']
                            fppos = hdulist[0].header['FPPOS']
                            if cenwave == 1096 and fppos == 1:
                                self.fix_cos_1096_dq(segment)
                        self.members.append(segment)

    def fix_cos_1096_dq(self, segment):
        print('Fixing COS 1096 FPPOS=1 data')
        segmentname = segment.data['segment']
        lastgoodwavelength = {'FUVA': 1239.75, 'FUVB': 1085.75}
        wavelength = segment.data['wavelength']
        good_indices = np.where(wavelength <= lastgoodwavelength[segmentname])
        lastgoodindex = good_indices[0][-1]
        segment.data['dq'][lastgoodindex+1:] |= 8
        return
        
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
        self.output_sumflux = np.zeros(self.nelements)
        self.output_sumweight = np.zeros(self.nelements)
        self.output_flux = np.zeros(self.nelements)
        self.output_errors = np.zeros(self.nelements)
        self.signal_to_noise = np.zeros(self.nelements)
        self.output_exptime = np.zeros(self.nelements)

        return wavegrid

    def wavelength_to_index(self, wavelength):
        index = (wavelength - self.min_wavelength) / self.delta_wavelength
        indices = [int(round(x)) for x in index]
        return indices

    def index_to_wavelength(self, index):
        wavelength = index * self.delta_wavelength + self.min_wavelength
        return wavelength

    def get_gross_counts(self, segment):
        pass

    def coadd(self):
        self.target = self.ull_targname()
        self.targ_ra, self.targ_dec = self.ull_coords()
        
        for segment in self.members:
            goodpixels = np.where((segment.data['dq'] & segment.sdqflags) == 0)
            wavelength = segment.data['wavelength'][goodpixels]
            indices = self.wavelength_to_index(wavelength)
            gross_counts = self.get_gross_counts(segment)
            weight = gross_counts[goodpixels]
            flux = segment.data['flux'][goodpixels]
            self.output_sumweight[indices] = self.output_sumweight[indices] + weight
            self.output_sumflux[indices] = self.output_sumflux[indices] + flux * weight
            self.output_exptime[indices] = self.output_exptime[indices] + segment.exptime
        nonzeros = np.where(self.output_sumweight != 0)
        if self.instrument == 'COS':
            # Using the variances (which only COS has) gives spikes in the error when the flux goes negative.
            self.output_sumweight[nonzeros] = np.where(self.output_sumweight[nonzeros] < 0.5, 0.5, self.output_sumweight[nonzeros])
        self.output_flux[nonzeros] = self.output_sumflux[nonzeros] / self.output_sumweight[nonzeros]
        # For the moment calculate errors from the gross counts
        self.output_errors[nonzeros] = np.sqrt(self.output_sumweight[nonzeros])
        self.signal_to_noise[nonzeros] = self.output_sumweight[nonzeros] / self.output_errors[nonzeros]
        self.output_errors[nonzeros] = np.abs(self.output_flux[nonzeros] / self.signal_to_noise[nonzeros])
        return

    def write(self, filename, overwrite=False, level="", version=""):
        
        # Table 1 - HLSP data
    
        # set up the header
        hdr1 = fits.Header()
        hdr1['EXTNAME'] = ('SCIENCE', 'Spectrum science arrays')
        hdr1['TIMESYS'] = ('UTC', 'Time system in use')
        hdr1['TIMEUNIT'] = ('s', 'Time unit for durations')
        hdr1['TREFPOS'] = ('GEOCENTER', 'Time reference position')

        mjd_beg = self.combine_keys("expstart", 1, "min")
        mjd_end = self.combine_keys("expend", 1, "max")
        dt_beg = Time(mjd_beg, format="mjd").datetime
        dt_end = Time(mjd_end, format="mjd").datetime
        hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S"), 'Date-time of first observation start')
        hdr1.add_blank('', after='TREFPOS')
        hdr1.add_blank('              / FITS TIME COORDINATE KEYWORDS', before='DATE-BEG')
    
        hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
        hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
        hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
        hdr1['XPOSURE'] = (self.combine_keys("exptime", 1, "sum"), '[s] Sum of exposure durations')
    
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
        hdr0['TELESCOP'] = (self.combine_keys("telescop", 0, "multi"), 'Telescope used to acquire data')
        hdr0['INSTRUME'] = (self.combine_keys("instrume", 0, "multi"), 'Instrument used to acquire data')
        hdr0.add_blank('', after='TELESCOP')
        hdr0.add_blank('              / SCIENCE INSTRUMENT CONFIGURATION', before='INSTRUME')
        hdr0['DETECTOR'] = (self.combine_keys("detector", 0, "multi"), 'Detector or channel used to acquire data')
        hdr0['DISPERSR'] = (self.combine_keys("opt_elem", 0, "multi"), 'Identifier of disperser')
        hdr0['CENWAVE'] = (self.combine_keys("cenwave", 0, "multi"), 'Central wavelength setting for disperser')
        hdr0['APERTURE'] = (self.combine_keys("aperture", 0, "multi"), 'Identifier of entrance aperture')
        hdr0['S_REGION'] = (self.obs_footprint(), 'Region footprint')
        hdr0['OBSMODE'] = (self.combine_keys("obsmode", 0, "multi"), 'Instrument operating mode (ACCUM | TIME-TAG)')
        hdr0['TARGNAME'] = self.targname[0]
        hdr0.add_blank(after='OBSMODE')
        hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

        hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
        hdr0['TARG_RA'] =  (self.targ_ra,  '[deg] Target right ascension')
        hdr0['TARG_DEC'] =  (self.targ_dec,  '[deg] Target declination')
        hdr0['PROPOSID'] = (self.combine_keys("proposid", 0, "multi"), 'Program identifier')
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
        hdr0['REFERENC'] = ('TBD', 'Bibliographic ID of primary paper')
    
        hdr0['CENTRWV'] = (self.combine_keys("centrwv", 0, "average"), 'Central wavelength of the data')
        hdr0.add_blank(after='REFERENC')
        hdr0.add_blank('           / ARCHIVE SEARCH KEYWORDS', before='CENTRWV')
        hdr0['MINWAVE'] = (self.combine_keys("minwave", 0, "min"), 'Minimum wavelength in spectrum')
        hdr0['MAXWAVE'] = (self.combine_keys("maxwave", 0, "max"), 'Maximum wavelength in spectrum')

        primary = fits.PrimaryHDU(header=hdr0)

        # Table 2 - individual product info
    
        # first set up header
        hdr2 = fits.Header()
        hdr2['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
        # set up the table columns
        cfn = fits.Column(name='FILENAME', array=np.array([h["filename"] for h in self.primary_headers]), format='A32')
        cpid = fits.Column(name='PROPOSID', array=np.array([h["proposid"] for h in self.primary_headers]), format='A32')
        ctel = fits.Column(name='TELESCOPE', array=np.array([h["telescop"] for h in self.primary_headers]), format='A32')
        cins = fits.Column(name='INSTRUMENT', array=np.array([h["instrume"] for h in self.primary_headers]), format='A32')
        cdet = fits.Column(name='DETECTOR', array=np.array([h["detector"] for h in self.primary_headers]), format='A32')
        cdis = fits.Column(name='DISPERSER', array=np.array([h["opt_elem"] for h in self.primary_headers]), format='A32')
        ccen = fits.Column(name='CENWAVE', array=np.array([h["cenwave"] for h in self.primary_headers]), format='A32')
        cap = fits.Column(name='APERTURE', array=np.array([h["aperture"] for h in self.primary_headers]), format='A32')
        csr = fits.Column(name='SPECRES', array=np.array([h["specres"] for h in self.primary_headers]), format='F8.1')
        ccv = fits.Column(name='CAL_VER', array=np.array([h["cal_ver"] for h in self.primary_headers]), format='A32')
        mjd_begs = np.array([h["expstart"] for h in self.first_headers])
        mjd_ends = np.array([h["expend"] for h in self.first_headers])
        mjd_mids = (mjd_ends + mjd_begs) / 2.
        cdb = fits.Column(name='MJD_BEG', array=mjd_begs, format='F15.9', unit='d')
        cdm = fits.Column(name='MJD_MID', array=mjd_mids, format='F15.9', unit='d')
        cde = fits.Column(name='MJD_END', array=mjd_ends, format='F15.9', unit='d')
        cexp = fits.Column(name='XPOSURE', array=np.array([h["exptime"] for h in self.first_headers]), format='F15.9', unit='s')
        cmin = fits.Column(name='MINWAVE', array=np.array([h["minwave"] for h in self.primary_headers]), format='F9.4', unit='Angstrom')
        cmax = fits.Column(name='MAXWAVE', array=np.array([h["maxwave"] for h in self.primary_headers]), format='F9.4', unit='Angstrom')
    
        cd2 = fits.ColDefs([cfn, cpid, ctel, cins, cdet, cdis, ccen, cap, csr, ccv, cdb, cdm, cde, cexp, cmin ,cmax])
    
        table2 = fits.BinTableHDU.from_columns(cd2, header=hdr2)
    
        # the HDUList:
        # 0 - empty data - 0th ext header
        # 1 - HLSP data - 1st ext header
        # 2 - individual product info - 2nd ext header
    
        hdul = fits.HDUList([primary, table1, table2])
    
        hdul.writeto(filename, overwrite=overwrite)

    def obs_footprint(self):
        radius = (2.5 / 2 / 3600)
        center_ra = self.targ_ra
        center_dec = self.targ_dec

        s_region = f"CIRCLE {center_ra} {center_dec} {radius}"
        return s_region

    def ull_targname(self):
        aliases = pd.read_json("pd_all_aliases.json", orient="split")
        targ_set = list(set([h["targname"] for h in self.primary_headers]))
        if len(targ_set) == 1:
            ull_targname = targ_set[0]
        else:
            ull_targname = self.primary_headers[0]["targname"]
        targ_matched = False
        for targ in self.targname:
            mask = aliases.apply(lambda row: row.astype(str).str.fullmatch(re.escape(targ)).any(), axis=1)
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
        
        master_list = pd.read_json("pd_targetinfo.json", orient="split")
        coords = master_list.loc[master_list["mast_targname"] == self.target][["ra", "dec"]].values
        if len(coords) != 0:
            return coords[0][0], coords[0][1]
        else:
            return avg_ra, avg_dec
                                      
    def combine_keys(self, key, hdrno, method):
        # Allowable methods are min, max, average, sum, multi
        if hdrno == 0:
            hdrs = self.primary_headers
        else:
            hdrs = self.first_headers

        vals = [h[key] for h in hdrs]
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


class STISSegmentList(SegmentList):

    def get_gross_counts(self, segment):
       exptime = segment.exptime
       gross = segment.data['gross']
       return np.abs(gross*exptime)


class COSSegmentList(SegmentList):

    def get_gross_counts(self, segment):
        try:
            gross = segment.data['variance_counts'] + segment.data['variance_bkg'] + segment.data['variance_flat']
            return gross
        except KeyError:
            gross = segment.data['gcounts']
            return gross


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
            transition_index_short = product_short.nelements
            transition_index_long = 0
        output_grating = product_short.grating + '-' + product_long.grating
        product_abutted = SegmentList(output_grating)
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
        if product_short.instrument == product_long.instrument:
            product_abutted.instrument = product_short.instrument
        else:
            product_abutted.instrument = product_short.instrument + '-' + product_long.instrument
        target_matched = False
        if product_short.target == product_long.target:
            product_abutted.target = product_short.target
            target_matched = True
            product_abutted.targname = [product_short.target, product_long.target]
        else:
            for target_name in product_short.targname:
                if target_name in product_long.targname:
                    product_abutted.target = product_abutted.ull_targname()
                    target_matched = True
                    product_abutted.targname = [target_name]
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
