import re
import glob
import os
import astropy
import pandas as pd
from astropy.io import fits
import numpy as np
import datetime
from datetime import datetime as dt
from astropy.time import Time

from ullyses_utils.parse_csv import parse_aliases
from ullyses_utils.ullyses_config import VERSION, CAL_VER

CODEDIR = os.path.dirname(__file__)
RED = "\033[1;31m"
RESET = "\033[0;0m"

class Ullyses():
    def __init__(self, files, hlspname, targname, ra, dec, level,
                 cal_ver=CAL_VER, version=VERSION, hlsp_type="spectral", 
                 overwrite=True, photfile=None):
        self.hlsp_type = hlsp_type
        if hlsp_type == "lcogt":
            assert photfile is not None, "Photometry file must be supplied for LCOGT data"
            self.photfile = photfile
            df = pd.read_csv(self.photfile, 
                    names=["filename", "mjdstart", "mjdend", "wl", "flux", "err"],
                    skiprows=[0], delim_whitespace=True)
            df = df.sort_values("mjdstart")
            self.photdf = df
        self.files = files
        self.primary_headers = []
        self.first_headers = []
        self.second_headers = []
        for item in self.files:
            self.primary_headers.append(fits.getheader(item))
            self.first_headers.append(fits.getheader(item, 1))
            if hlsp_type == "drizzled":
                self.second_headers.append(fits.getheader(item, 2))
        self.targname = targname
        self.targ_ra = ra
        self.targ_dec = dec
        self.hlspname = hlspname
        self.cal_ver = cal_ver
        self.version = version
        self.level = level
        self.overwrite = overwrite


    def make_hdrs_and_prov(self):
        if self.hlsp_type == "spectral":
            self.make_spectral_hdr0()
            self.make_spectral_hdr1()
            self.make_spectral_prov_ext()
        elif self.hlsp_type == "imaging":
            self.make_imaging_hdr0()
            self.make_imaging_hdr1()
            self.make_imaging_prov_ext()
        elif self.hlsp_type == "lcogt":
            self.make_lcogt_timeseries_hdr0()
            self.make_lcogt_timeseries_data_ext()
            self.make_lcogt_timeseries_prov_ext()
        elif self.hlsp_type == "drizzled":
            assert len(self.files) == 1, f"{len(self.files)} provided, can only handle 1"
            self.make_imaging_hdr0()
            self.make_drizzled_data_ext()
            self.make_drizzled_wgt_ext()
            self.make_drizzled_prov_ext()
        else:
            print(f"ERROR: HLSP type not '{self.hlsp_type}' recognized. Must be 'spectral', 'imaging', 'drizzled', or 'lcogt'")

        
    def write_file(self):
        if self.hlsp_type == "spectral":
            pass
        elif self.hlsp_type == "imaging":
            pass
        elif self.hlsp_type == "lcogt":
            primary = fits.PrimaryHDU(header=self.hdr0)
            ext1 = self.hdu1
            prov = self.prov_hdu
            hdulist = fits.HDUList([primary, ext1, prov])
        elif self.hlsp_type == "drizzled":
            assert len(self.files) == 1, f"{len(self.files)} provided, can only handle 1"
            primary = fits.PrimaryHDU(header=self.hdr0)
            ext1 =  self.hdu1
            ext2 = self.hdu2
            prov = self.prov_hdu
            hdulist = fits.HDUList([primary, ext1, ext2, prov])
        else:
            print(f"ERROR: HLSP type not '{self.hlsp_type}' recognized. Must be 'spectral', 'imaging', 'drizzled', or 'lcogt'")
        hdulist.writeto(self.hlspname, overwrite=self.overwrite)
        print(f"Wrote {self.hlspname}")


    def make_spectral_hdr0(self):
        hdr0 = fits.Header()
        hdr0['EXTEND'] = ('T', 'FITS file may contain extensions')
        hdr0['NEXTEND'] = 3
        hdr0['FITS_VER'] = 'Definition of the Flexible Image Transport System (FITS) v4.0 https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf'
        hdr0['FITS_SW'] = ('astropy.io.fits v' + astropy.__version__, 'FITS file creation software')
        hdr0['ORIGIN'] = ('Space Telescope Science Institute', 'FITS file originator')
        hdr0['DATE'] = (str(datetime.date.today()), 'Date this file was written')
        hdr0['FILENAME'] = (os.path.basename(self.hlspname), 'Name of this file')
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
        hdr0['TARGNAME'] = self.targname
        hdr0.add_blank(after='OBSMODE')
        hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

        hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
        hdr0['TARG_RA'] =  (self.targ_ra,  '[deg] Target right ascension')
        hdr0['TARG_DEC'] =  (self.targ_dec,  '[deg] Target declination')
        hdr0['PROPOSID'] = (self.combine_keys("proposid", "multi"), 'Program identifier')
        hdr0.add_blank(after='TARG_DEC')
        hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
        hdr0['CAL_VER'] = (f'ULLYSES Cal {self.cal_ver}', 'HLSP processing software version')
        hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
        hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                        'Name ID of this HLSP collection')
        hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead')
        hdr0['HLSP_VER'] = (self.version, 'HLSP data release version identifier')
        hdr0['HLSP_LVL'] = (self.level, 'ULLYSES HLSP Level')
        hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
        hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
        hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')  
                                                                                                                           
        hdr0['CENTRWV'] = (self.combine_keys("centrwv", "average"), 'Central wavelength of the data')                      
        hdr0.add_blank(after='REFERENC')                                                                                   
        hdr0.add_blank('           / ARCHIVE SEARCH KEYWORDS', before='CENTRWV')                                           
        hdr0['MINWAVE'] = (self.combine_keys("minwave", "min"), 'Minimum wavelength in spectrum')                          
        hdr0['MAXWAVE'] = (self.combine_keys("maxwave", "max"), 'Maximum wavelength in spectrum')
    
        self.hdr0 = hdr0 


    def make_imaging_hdr0(self):
        hdr0 = fits.Header()
        hdr0['EXTEND'] = ('T', 'FITS file may contain extensions')
        hdr0['NEXTEND'] = 3
        hdr0['FITS_VER'] = 'Definition of the Flexible Image Transport System (FITS) v4.0 https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf'
        hdr0['FITS_SW'] = ('astropy.io.fits v' + astropy.__version__, 'FITS file creation software')
        hdr0['ORIGIN'] = ('Space Telescope Science Institute', 'FITS file originator')
        hdr0['DATE'] = (str(datetime.date.today()), 'Date this file was written')
        hdr0['FILENAME'] = (os.path.basename(self.hlspname), 'Name of this file')
        telescop = self.combine_keys("telescop", "multi")
        hdr0['TELESCOP'] = (telescop, 'Telescope used to acquire data')
        instrume = self.combine_keys("instrume", "multi")
        hdr0['INSTRUME'] = (instrume, 'Instrument used to acquire data')
        hdr0.add_blank('', after='TELESCOP')
        hdr0.add_blank('              / SCIENCE INSTRUMENT CONFIGURATION', before='INSTRUME')
        hdr0['DETECTOR'] = (self.combine_keys("detector", "multi"), 'Detector or channel used to acquire data')
        hdr0['FILTER'] = (self.combine_keys("filter", "multi"), 'Identifier of disperser')
        hdr0['APERTURE'] = (self.combine_keys("aperture", "multi"), 'Identifier of entrance aperture')
        if telescop == "HST":
            hdr0["FGSLOCK"] = (self.combine_keys("fgslock", "multi"), "Commanded FGS lock (FINE,COARSE,GYROS,UNKNOWN)")
            hdr0['GYROMODE'] = (self.combine_keys("gyromode", "multi"), "Number of gyros scheduled, T=3+OBAD")
            if instrume == "WFC3":
                hdr0["FLASHDUR"] = (self.combine_keys("flashdur", "multi"), "Post flash exposure time in seconds: 0.1 to 409.5")
                hdr0["FLASHCUR"] = (self.combine_keys("flashcur", "multi"), "Post flash current (zero, low, medium, high)")
                hdr0["FLASHLVL"] = (self.combine_keys("flashlvl", "multi"), "Post flash requested flash level")
                hdr0["FLASHSTA"] = (self.combine_keys("flashsta", "multi"), "Post flash status: SUCCESSFUL, ABORTED, NOT PERFORMED")
        hdr0['OBSMODE'] = (self.combine_keys("obsmode", "multi"), 'Instrument operating mode (ACCUM | TIME-TAG)')
        hdr0['TARGNAME'] = self.targname
        hdr0.add_blank(after='OBSMODE')
        hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

        hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
        hdr0['TARG_RA'] =  (self.targ_ra,  '[deg] Target right ascension')
        hdr0['TARG_DEC'] =  (self.targ_dec,  '[deg] Target declination')
        hdr0['PROPOSID'] = (self.combine_keys("proposid", "multi"), 'Program identifier')
        hdr0.add_blank(after='TARG_DEC')
        hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
        hdr0['CAL_VER'] = (f'ULLYSES Cal {self.cal_ver}', 'HLSP processing software version')
        hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
        hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                        'Name ID of this HLSP collection')
        hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead')
        hdr0['HLSP_VER'] = (self.version, 'HLSP data release version identifier')
        hdr0['HLSP_LVL'] = (self.level, 'ULLYSES HLSP Level')
        hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
        hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
        hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')  
    
        self.hdr0 = hdr0
    

    def make_lcogt_timeseries_hdr0(self):
        hdr0 = fits.Header()
        hdr0['EXTEND'] = ('T', 'FITS file may contain extensions')
        hdr0['NEXTEND'] = 3
        hdr0['FITS_VER'] = 'Definition of the Flexible Image Transport System (FITS) v4.0 https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf'
        hdr0['FITS_SW'] = ('astropy.io.fits v' + astropy.__version__, 'FITS file creation software')
        hdr0['ORIGIN'] = ('Space Telescope Science Institute', 'FITS file originator')
        hdr0['DATE'] = (str(datetime.date.today()), 'Date this file was written')
        hdr0['FILENAME'] = (os.path.basename(self.hlspname), 'Name of this file')
        telescop = self.combine_keys("telescop", "multi", constant="LCOGT")
        hdr0['TELESCOP'] = (telescop, 'Telescope used to acquire data')
        instrume = self.combine_keys("instrume", "multi", dict_key="LCOGT")
        hdr0['INSTRUME'] = (instrume, 'Instrument used to acquire data')
        hdr0.add_blank('', after='TELESCOP')
        hdr0.add_blank('              / SCIENCE INSTRUMENT CONFIGURATION', before='INSTRUME')
        hdr0['DETECTOR'] = (self.combine_keys("detector", "multi", constant="LCOGT"), 'Detector or channel used to acquire data')
        hdr0['FILTER'] = (self.combine_keys("filter", "multi", dict_key="LCOGT"), 'Identifier of disperser')
        hdr0['APERTURE'] = (self.combine_keys("aperture", "multi", constant="LCOGT"), 'Identifier of entrance aperture')
        hdr0['OBSMODE'] = (self.combine_keys("obsmode", "multi", constant="IMAGING"), 'Instrument operating mode (ACCUM | TIME-TAG)')
        hdr0['TARGNAME'] = self.targname
        hdr0.add_blank(after='OBSMODE')
        hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

        hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
        hdr0['TARG_RA'] =  (self.targ_ra,  '[deg] Target right ascension')
        hdr0['TARG_DEC'] =  (self.targ_dec,  '[deg] Target declination')
        hdr0['PROPOSID'] = (self.combine_keys("proposid", "multi", dict_key="LCOGT"), 'Program identifier')
        hdr0.add_blank(after='TARG_DEC')
        hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
        hdr0['CAL_VER'] = (f'ULLYSES Cal {self.cal_ver}', 'HLSP processing software version')
        hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
        hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                        'Name ID of this HLSP collection')
        hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead')
        hdr0['HLSP_VER'] = (self.version, 'HLSP data release version identifier')
        hdr0['HLSP_LVL'] = (self.level, 'ULLYSES HLSP Level')
        hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
        hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
        hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')  
        # CENTRWV?

        self.hdr0 = hdr0


    def make_spectral_hdr1(self):
        hdr1 = fits.Header()
        hdr1['EXTNAME'] = ('SCIENCE', 'Spectrum science arrays')
        hdr1['TIMESYS'] = ('UTC', 'Time system in use')
        hdr1['TIMEUNIT'] = ('s', 'Time unit for durations')
        hdr1['TREFPOS'] = ('GEOCENTER', 'Time reference position')

        mjd_beg = self.combine_keys("expstart", "min")
        mjd_end = self.combine_keys("expend", "max")
        dt_beg = Time(mjd_beg, format="mjd").datetime
        dt_end = Time(mjd_end, format="mjd").datetime
        hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S%.f"), 'Date-time of first observation start')
        hdr1.add_blank('', after='TREFPOS')
        hdr1.add_blank('              / FITS TIME COORDINATE KEYWORDS', before='DATE-BEG')

        hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
        hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
        hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
        hdr1['XPOSURE'] = (self.combine_keys("exptime", "sum"), '[s] Sum of exposure durations')
    
        self.hdr1 = hdr1
                                                                                                                           
    def make_imaging_hdr1(self):
        hdr1 = fits.Header()
        hdr1['EXTNAME'] = ('SCIENCE', 'Image array')
        hdr1['TIMESYS'] = ('UTC', 'Time system in use')
        hdr1['TIMEUNIT'] = ('s', 'Time unit for durations')
        hdr1['TREFPOS'] = ('GEOCENTER', 'Time reference position')

        mjd_beg = self.combine_keys("expstart", "min", "WFC3")
        mjd_end = self.combine_keys("expend", "max", "WFC3")
        dt_beg = Time(mjd_beg, format="mjd").datetime
        dt_end = Time(mjd_end, format="mjd").datetime
        hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S"), 'Date-time of first observation start')
        hdr1.add_blank('', after='TREFPOS')
        hdr1.add_blank('              / FITS TIME COORDINATE KEYWORDS', before='DATE-BEG')

        hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
        hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
        hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
        hdr1['XPOSURE'] = (self.combine_keys("exptime", "sum", "WFC3"), '[s] Sum of exposure durations')
    
        self.hdr1 = hdr1
    

    def make_lcogt_timeseries_data_ext(self):
        hdr1 = fits.Header()
        hdr1['EXTNAME'] = ('SCIENCE', 'Spectrum science arrays')
        hdr1['TIMESYS'] = ('UTC', 'Time system in use')
        hdr1['TIMEUNIT'] = ('s', 'Time unit for durations')
        hdr1['TREFPOS'] = ('GEOCENTER', 'Time reference position')

        mjd_beg = self.photdf.iloc[0]["mjdstart"]
        mjd_end = self.photdf.iloc[-1]["mjdend"]
        dt_beg = Time(mjd_beg, format="mjd").datetime
        dt_end = Time(mjd_end, format="mjd").datetime
        hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S"), 'Date-time of first observation start')
        hdr1.add_blank('', after='TREFPOS')
        hdr1.add_blank('              / FITS TIME COORDINATE KEYWORDS', before='DATE-BEG')

        hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
        hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
        hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
        hdr1['XPOSURE'] = (self.combine_keys("exptime", "sum", "LCOGT"), '[s] Sum of exposure durations')
    
        self.hdr1 = hdr1

        nrows = len(self.photdf["mjdstart"])
        ncols = len(set(self.photdf["wl"]))
        output_flux = np.zeros((nrows, ncols))
        output_err = np.zeros((nrows, ncols))
        output_time0 = self.photdf["mjdstart"].values
        output_time1 = self.photdf["mjdend"].values
        output_wl = np.array(list(set(self.photdf["wl"].values)))
        output_wl.sort()
        for i in range(len(self.photdf)):
            wlind = np.where(output_wl == self.photdf.iloc[i]["wl"])[0]
            output_flux[i, wlind] = self.photdf.iloc[i]["flux"]
            output_err[i, wlind] = self.photdf.iloc[i]["err"]

        npixels = nrows * ncols
        array_dimensions = '(' + str(ncols) + ', ' + str(nrows) + ')'
        columns = []
        columns.append(fits.Column(name='MJDSTART', format=str(nrows)+'D', unit='Day'))
        columns.append(fits.Column(name='MJDEND', format=str(nrows)+'D', unit='Day'))
        columns.append(fits.Column(name='WAVELENGTH', format=str(ncols)+'E', 
                                   unit='Angstrom'))
        columns.append(fits.Column(name='FLUX', format=str(npixels)+'E', 
                                   dim=array_dimensions, unit='erg /s /cm**2 /Angstrom'))
        columns.append(fits.Column(name='ERROR', format=str(npixels)+'E', 
                                   dim=array_dimensions, unit='erg /s /cm**2 /Angstrom'))
        cd = fits.ColDefs(columns)
        table1 = fits.BinTableHDU.from_columns(cd, header=hdr1, nrows=1)
        table1.data[0]["wavelength"] = output_wl
        table1.data[0]["mjdstart"] = output_time0
        table1.data[0]["mjdend"] = output_time1
        table1.data[0]["flux"] = output_flux
        table1.data[0]["error"] = output_err

        self.data1 = cd
        self.hdu1 = table1


    def make_drizzled_data_ext(self):
        hdr1 = self.first_headers[0]
        for key in ["ROOTNAME", "EXPNAME", "HDRNAME", "IDCTAB", "FITNAMEB"]:
            try:
                del hdr1[key]
            except KeyError:
                pass
        hdr1["DRIZSCAL"] = (self.primary_headers[0]["D001SCAL"], "Drizzle: pixel size (arcsec) of output image")
        hdr1["DRIZISCL"] = (self.primary_headers[0]["D001ISCL"], "Drizzle: default IDCTAB pixel size(arcsec)")
        hdr1["DRIZPIXF"] = (self.primary_headers[0]["D001PIXF"], "Drizzle: linear size of drop")
        mjd_beg = self.combine_keys("expstart", "min", "WFC3")
        mjd_end = self.combine_keys("expend", "max", "WFC3")
        dt_beg = Time(mjd_beg, format="mjd").datetime
        dt_end = Time(mjd_end, format="mjd").datetime
        hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S"), 'Date-time of first observation start')

        hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
        hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
        hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
        hdr1['XPOSURE'] = (self.combine_keys("exptime", "sum", "WFC3"), '[s] Sum of exposure durations')
        
        data1 = fits.getdata(self.files[0], 1)

        self.hdr1 = hdr1
        self.data1 = data1
        self.hdu1 = fits.ImageHDU(data1, header=hdr1) 
    
    
    def make_drizzled_wgt_ext(self):
        hdr2 = self.second_headers[0]
        del hdr2["ROOTNAME"]
        del hdr2["EXPNAME"]
        hdr2["DRIZSCAL"] = (self.primary_headers[0]["D001SCAL"], "Drizzle: pixel size (arcsec) of output image")
        hdr2["DRIZISCL"] = (self.primary_headers[0]["D001ISCL"], "Drizzle: default IDCTAB pixel size(arcsec)")
        hdr2["DRIZPIXF"] = (self.primary_headers[0]["D001PIXF"], "Drizzle: linear size of drop")
        
        data2 = fits.getdata(self.files[0], 2)

        self.hdr2 = hdr2
        self.data2 = data2
        self.hdu2 = fits.ImageHDU(data2, header=hdr2) 
          

    def make_spectral_prov_ext(self):
        hdr = fits.Header()
        hdr['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
        # set up the table columns
        cfn = fits.Column(name='FILENAME', array=self.combine_keys("filename", "arr"), 
            format='A40')
        cpid = fits.Column(name='PROPOSID', array=self.combine_keys("proposid", "arr"), 
            format='A32')
        ctel = fits.Column(name='TELESCOPE', array=self.combine_keys("telescop", "arr"), 
            format='A32')
        cins = fits.Column(name='INSTRUMENT', array=self.combine_keys("instrume", "arr"), 
            format='A32')
        cdet = fits.Column(name='DETECTOR', array=self.combine_keys("detector", "arr"), 
            format='A32')
        cdis = fits.Column(name='DISPERSER', array=self.combine_keys("opt_elem", "arr"), 
            format='A32')
        ccen = fits.Column(name='CENWAVE', array=self.combine_keys("cenwave", "arr"), 
            format='A32')
        cap = fits.Column(name='APERTURE', array=self.combine_keys("aperture", "arr"), 
            format='A32')
        csr = fits.Column(name='SPECRES', array=self.combine_keys("specres", "arr"), 
            format='F8.1')
        ccv = fits.Column(name='CAL_VER', array=self.combine_keys("cal_ver", "arr"), 
            format='A32')
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
        table2 = fits.BinTableHDU.from_columns(cd2, header=hdr)

        self.prov_hdr = hdr
        self.prov_data = cd2
        self.prov_hdu = table2
    

    def make_imaging_prov_ext(self):
        hdr = fits.Header()
        hdr['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
        # set up the table columns
        cfn = fits.Column(name='FILENAME', array=self.combine_keys("filename", "arr"), 
            format='A40')
        cpid = fits.Column(name='PROPOSID', array=self.combine_keys("proposid", "arr"), 
            format='A32')
        ctel = fits.Column(name='TELESCOPE', array=self.combine_keys("telescop", "arr"), 
            format='A32')
        cins = fits.Column(name='INSTRUMENT', array=self.combine_keys("instrume", "arr"), 
            format='A32')
        cdet = fits.Column(name='DETECTOR', array=self.combine_keys("detector", "arr"), 
            format='A32')
        cfil = fits.Column(name='FILTER', array=self.combine_keys("filter", "arr"), 
            format='A32')
        cap = fits.Column(name='APERTURE', array=self.combine_keys("aperture", "arr"), 
            format='A32')
        ccv = fits.Column(name='CAL_VER', array=self.combine_keys("cal_ver", "arr"), 
            format='A32')
        mjd_begs = self.combine_keys("expstart", "arr", "WFC3")
        mjd_ends = self.combine_keys("expend", "arr", "WFC3")
        mjd_mids = (mjd_ends + mjd_begs) / 2.
        cdb = fits.Column(name='MJD_BEG', array=mjd_begs, format='F15.9', unit='d')
        cdm = fits.Column(name='MJD_MID', array=mjd_mids, format='F15.9', unit='d')
        cde = fits.Column(name='MJD_END', array=mjd_ends, format='F15.9', unit='d')
        cexp = fits.Column(name='XPOSURE', array=self.combine_keys("exptime", "arr", "WFC3"), format='F15.9', unit='s')

        cd2 = fits.ColDefs([cfn, cpid, ctel, cins, cdet, cfil, cap, ccv, cdb, cdm, cde, cexp])
        table2 = fits.BinTableHDU.from_columns(cd2, header=hdr)

        self.prov_hdr = hdr
        self.prov_data = cd2
        self.prov_hdu = table2
    

    def make_lcogt_timeseries_prov_ext(self):
        hdr = fits.Header()
        hdr['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
        # set up the table columns
        files = [os.path.basename(x) for x in self.files]
        cfn = fits.Column(name='FILENAME', array=files,
            format='A40')
        cpid = fits.Column(name='PROPOSID', array=self.combine_keys("proposid", "arr", "LCOGT"), 
            format='A32')
        ctel = fits.Column(name='TELESCOPE', array=self.combine_keys("telescop", "arr", "LCOGT"), 
            format='A32')
        cins = fits.Column(name='INSTRUMENT', array=self.combine_keys("instrume", "arr", "LCOGT"), 
            format='A32')
        cdet = fits.Column(name='DETECTOR', array=self.combine_keys("detector", "arr", constant="LCOGT"), 
            format='A32')
        cfil = fits.Column(name='FILTER', array=self.combine_keys("filter", "arr", "LCOGT"), 
            format='A32')
        cap = fits.Column(name='APERTURE', array=self.combine_keys("aperture", "arr", constant="LCOGT"), 
            format='A32')
        ccv = fits.Column(name='CAL_VER', array=self.combine_keys("cal_ver", "arr", "LCOGT"), 
            format='A32')
        
        mjd_begs = self.photdf["mjdstart"].to_numpy()
        mjd_ends = self.photdf["mjdend"].to_numpy()
        mjd_mids = (mjd_ends + mjd_begs) / 2.
        cdb = fits.Column(name='MJD_BEG', array=mjd_begs, format='F15.9', unit='d')
        cdm = fits.Column(name='MJD_MID', array=mjd_mids, format='F15.9', unit='d')
        cde = fits.Column(name='MJD_END', array=mjd_ends, format='F15.9', unit='d')
        cexp = fits.Column(name='XPOSURE', array=self.combine_keys("exptime", "arr", "LCOGT"), format='F15.9', unit='s')

        cd2 = fits.ColDefs([cfn, cpid, ctel, cins, cdet, cfil, cap, ccv, cdb, cdm, cde, cexp])
        table2 = fits.BinTableHDU.from_columns(cd2, header=hdr)

        self.prov_hdr = hdr
        self.prov_data = cd2
        self.prov_hdu = table2


    def make_drizzled_prov_ext(self):
        hdr = fits.Header()
        data = fits.getdata(self.files[0], 4)
        inds = []
        files = []
        for i,f in enumerate(data["filename"]):
            if f not in files:
                files.append(f)
                inds.append(i)

        hdr['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
        # set up the table columns
        cfn = fits.Column(name='FILENAME', array=data["filename"][inds], 
            format='A40')
        cpid = fits.Column(name='PROPOSID', array=data["proposid"][inds],
            format='A32')
        ctel = fits.Column(name='TELESCOPE', array=data["telescop"][inds],
            format='A32')
        cins = fits.Column(name='INSTRUMENT', array=data["instrume"][inds],
            format='A32')
        cdet = fits.Column(name='DETECTOR', array=data["detector"][inds],
            format='A32')
        cfil = fits.Column(name='FILTER', array=data["filter"][inds],
            format='A32')
        cap = fits.Column(name='APERTURE', array=data["aperture"][inds],
            format='A32')
        ccv = fits.Column(name='CAL_VER', array=data["cal_ver"][inds],
            format='A32')
        mjd_begs = data["expstart"][inds]
        mjd_ends = data["expend"][inds] 
        mjd_mids = (mjd_ends + mjd_begs) / 2.
        cdb = fits.Column(name='MJD_BEG', array=mjd_begs, format='F15.9', unit='d')
        cdm = fits.Column(name='MJD_MID', array=mjd_mids, format='F15.9', unit='d')
        cde = fits.Column(name='MJD_END', array=mjd_ends, format='F15.9', unit='d')
        cexp = fits.Column(name='XPOSURE', array=data["exptime"][inds], format='F15.9', unit='s')

        cd2 = fits.ColDefs([cfn, cpid, ctel, cins, cdet, cfil, cap, ccv, cdb, cdm, cde, cexp])
        table2 = fits.BinTableHDU.from_columns(cd2, header=hdr)
        
        self.prov_hdr = hdr
        self.prov_data = cd2
        self.prov_hdu = table2


    def combine_keys(self, key, method, dict_key=None, constant=None):
        keymap= {"HST": {"expstart": ("expstart", 1),
                         "expend": ("expend", 1),
                         "exptime": ("exptime", 1),
                         "telescop": ("telescop", 0),
                         "instrume": ("instrume", 0),
                         "detector": ("detector", 0),
                         "opt_elem": ("opt_elem", 0),
                         "filter": ("filter", 0),
                         "fgslock": ("fgslock", 0),
                         "gyromode": ("gyromode", 0),
                         "flashdur": ("flashdur", 0),
                         "flashcur": ("flashcur", 0),
                         "flashlvl": ("flashlvl", 0),
                         "flashsta": ("flashsta", 0),
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
                "WFC3": {"expstart": ("expstart", 0),
                         "expend": ("expend", 0),
                         "exptime": ("exptime", 0)},
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
                         "cal_ver": ("cf_vers", 0)},
               "LCOGT": {"expstart": ("date-obs", 1),
                         "expend": ("exptime", 1),
                         "exptime": ("exptime", 1),
                         "telescop": ("telescop", 1),
                         "instrume": ("instrume", 1),
                         "detector": ("telescop", 1),
                         "opt_elem": ("telescop", 1),
                         "proposid": ("propid", 1),
                         "filename": ("origname", 1),
                         "filter": ("filter", 1),
                         "cal_ver": ("pipever", 1)}
                         }

        if constant is not None:
            vals = [constant for x in self.primary_headers]
        else:
            vals = []
            for i in range(len(self.primary_headers)):
                if dict_key is None:
                    tel = self.primary_headers[i]["telescop"]
                else:
                    tel = dict_key
                actual_key = keymap[tel][key][0]
                hdrno = keymap[tel][key][1]
                if hdrno == 0:
                    val = self.primary_headers[i][actual_key]
                else:
                    val = self.first_headers[i][actual_key]
                if tel == "FUSE" and key == "filename":
                    val = val.replace(".fit", "_vo.fits")
                if tel == "LCOGT":
                    dto = dt.strptime(self.first_headers[i]["date-obs"], "%Y-%m-%dT%H:%M:%S.%f")
                    t = Time(dto, format="datetime")
                    mjdstart = t.mjd
                    if key == "expstart":
                        val = mjdstart
                    if key == "expend":
                        exptime = self.first_headers[i]["exptime"]
                        val = mjdstart + (exptime/86400) # seconds in a day
                    elif key == "telescop":
                        val = f"LCOGT-{val}"

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
        radius = (2.5 / 2 / 3600)
        center_ra = self.targ_ra
        center_dec = self.targ_dec


    def filter_to_wl(self, filts):
        conv = {"V": 5500, "ip": 7718.28}
        wl = np.array([conv[x] for x in filts])
        return wl


    def mag_to_flux(self, mags, filts):
        """
        Vega magnitude V zero points from Bessell et al. (1998) 
        AB magnitude i' zero points from http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=SLOAN&asttype=
        """
        zpts = {"V": 363.1e-11, "ip": 1.27e-9} # Flambda zero points
        filt_zpts = np.array([zpts[x] for x in filts])
        flux = (10 ** (-mags / 2.5)) * filt_zpts
        return flux

