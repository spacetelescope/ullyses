import pandas as pd
from collections import defaultdict
import argparse
import os
import glob
import re
import datetime
from datetime import datetime as dt
import numpy as np

import astropy
from astropy.io import fits
from astropy.time import Time

from ullyses.coadd import COSSegmentList, STISSegmentList, FUSESegmentList, CCDSegmentList
from ullyses.coadd import abut, SegmentList
from ullyses.combine_header_keys import KeyBlender
import ullyses_utils
from . import __release__, __version__
from ullyses_utils import parse_csv, match_aliases
from hasp.grating_priority import create_level4_products

RED = "\033[1;31m"
RESET = "\033[0;0m"

ULLYSES_GRATING_PRIORITIES = {'COS/G130M': {'minwave': 900, 'maxwave': 1470, 'priority': 1},
                              'FUSE/FUSE': {'minwave': 912, 'maxwave': 1179.9, 'priority': 2},
                              'STIS/E140M': {'minwave': 1143.0, 'maxwave': 1727.2, 'priority': 3},
                              'COS/G160M': {'minwave': 1342, 'maxwave': 1800, 'priority': 4},
                              'STIS/E140H': {'minwave': 1141.1, 'maxwave': 1687.9, 'priority': 5},
                              'STIS/G140M': {'minwave': 1145.1, 'maxwave': 1741.9, 'priority': 6},
                              'STIS/E230M': {'minwave': 1606.7, 'maxwave': 3119.2, 'priority': 7},
                              'STIS/E230H': {'minwave': 1629.0, 'maxwave': 3156.0, 'priority': 8},
                              'STIS/G230M': {'minwave': 1641.8, 'maxwave': 3098.2, 'priority': 9},
                              'COS/G140L': {'minwave': 901, 'maxwave': 2150, 'priority': 10},
                              'STIS/G230MB': {'minwave': 1635.0, 'maxwave': 3184.5, 'priority': 11},
                              'COS/G185M': {'minwave': 1664, 'maxwave': 2134, 'priority': 12},
                              'COS/G225M': {'minwave': 2069, 'maxwave': 2526, 'priority': 13},
                              'COS/G285M': {'minwave': 2474, 'maxwave': 3221, 'priority': 14},
                              'STIS/G140L': {'minwave': 1138.4, 'maxwave': 1716.4, 'priority': 15},
                              'STIS/G430M': {'minwave': 3021.9, 'maxwave': 5610.1, 'priority': 16},
                              'STIS/G230L': {'minwave': 1582.0, 'maxwave': 3158.7, 'priority': 17},
                              'STIS/G230LB': {'minwave': 1667.1, 'maxwave': 3071.6, 'priority': 18},
                              'COS/G230L': {'minwave': 1650, 'maxwave': 3200, 'priority': 19},
                              'STIS/G750M': {'minwave': 5464.6, 'maxwave': 10645.1, 'priority': 20},
                              'STIS/G430L': {'minwave': 2895.9, 'maxwave': 5704.4, 'priority': 21},
                              'STIS/G750L': {'minwave': 5261.3, 'maxwave': 10252.3, 'priority': 22},
                      }

'''
This wrapper goes through each target folder in the ullyses data directory and find
the data and which gratings are present. This info is then fed into coadd.py.
'''

class Ullyses_SegmentList(KeyBlender, SegmentList):
    """This class is a mixin to add the project-specific write function and functions
    to populate the target name and coordinates

    """

    # Segments that should be ignored for each cenwave.
    # For ULLYSES this is G230L/2635/NUVC and G230L/2950/NUVC
#    self.bad_segments = {2635: 'NUVC', 2950: 'NUVC'}
    SegmentList.bad_segments = {2635: 'NUVC', 2950: 'NUVC'}


    def write(self, filename, overwrite=False, level="", version=""):

        # If the target is a ULLYSES target, use the official
        # target name and coords
        self.target = self.get_targname()
        self.targ_ra, self.targ_dec, self.coord_epoch = self.get_coords()

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
        all_comments = self.combine_keys("comment", "comment")
        for comment in all_comments:
            hdr1['COMMENT'] = (comment, "Calibration and/or quality comment")

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
        hdr0['TARGNAME'] = self.get_targname("target_name_ullyses")
        hdr0.add_blank(after='OBSMODE')
        hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

        hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
        hdr0['G_EPOCH'] =  (self.coord_epoch,  'Epoch of GAIA coordinates')
        hdr0['TARG_RA'] =  (self.targ_ra,  '[deg] Target right ascension')
        hdr0['TARG_DEC'] =  (self.targ_dec,  '[deg] Target declination')
        hdr0['PROPOSID'] = (self.combine_keys("proposid", "multi"), 'Program identifier')
        hdr0.add_blank(after='TARG_DEC')
        hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
        hdr0['CAL_VER'] = (f'ULLYSES Cal {__version__}', 'HLSP processing software version')
        hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
        hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                        'Name ID of this HLSP collection')
        hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead')
        hdr0['HLSP_VER'] = (version,'HLSP data release version identifier')
        hdr0['HLSP_LVL'] = (level, 'ULLYSES HLSP Level')
        hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
        hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
        hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')

        minwave = self.combine_keys("minwave", "min")
        maxwave = self.combine_keys("maxwave", "max")
        centrwv = ((maxwave - minwave)/2.) + minwave
        hdr0['CENTRWV'] = (centrwv, 'Central wavelength of the data')

        hdr0.add_blank(after='REFERENC')
        hdr0.add_blank('           / ARCHIVE SEARCH KEYWORDS', before='CENTRWV')
        hdr0['MINWAVE'] = (minwave, 'Minimum wavelength in spectrum')
        hdr0['MAXWAVE'] = (maxwave, 'Maximum wavelength in spectrum')

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

        ## sort the provenance by MJD start time
        time_sort = np.argsort(cdb.array) # indices
        # sort each of the columns using these indices
        col_in_order = [cfn, cpid, ctel, cins, cdet, cdis, ccen, cap, csr, ccv,
                        cdb, cdm, cde, cexp, cmin, cmax]
        sorted_coldef = []
        for col in col_in_order:
            col.array = col.array[time_sort]
            sorted_coldef.append(col)

        # turn into a ColDef to feed into the HDU
        cd2 = fits.ColDefs(sorted_coldef)

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
        self.targ_ra, self.targ_dec, self.coord_epoch = self.get_coords()
        radius = (2.5 / 2 / 3600)
        center_ra = self.targ_ra
        center_dec = self.targ_dec

        s_region = f"CIRCLE {center_ra} {center_dec} {radius}"
        return s_region

    def get_targname(self, targcol="target_name_hlsp"):
        aliases_file = ullyses_utils.__path__[0] + '/data/target_metadata/ullyses_aliases.csv'
        aliases = pd.read_csv(aliases_file)
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
                ull_targname = aliases[mask][targcol].values[0]
                break
        if targ_matched is False:
            print(f"{RED}WARNING: Could not match target name {ull_targname} to ULLYSES alias list{RESET}")
        return ull_targname

    def get_coords(self):
        matched = False
        ras = list(set([h["ra_targ"] for h in self.primary_headers]))
        decs = list(set([h["dec_targ"] for h in self.primary_headers]))
        ra = np.average(ras)
        dec = np.average(decs)
        epoch = "UNKNOWN"
        if self.target == "":
            print(f"{RED}WARNING: Could not determine coordinates{RESET}")
            return ra, dec, epoch
        csvs, metadata_dfs = parse_csv.parse_database_csv("all")
        ull_alias = match_aliases.match_aliases(self.target, "target_name_ullyses")
        for df in metadata_dfs:
            df['target_name_ullyses'] = df['target_name_ullyses'].str.upper()
            row = df.loc[df.target_name_ullyses == ull_alias]
            if len(row) == 1:
                matched = True
                ra = row.targ_ra.values[0]
                dec = row.targ_dec.values[0]
                epoch = float(row.coordinate_epoch.values[0])
                break
        if matched is False:
            print(f"{RED}WARNING: Could not determine coordinates for {ull_alias}{RESET}")
        return ra, dec, epoch

    def add_hasp_attributes(self):
        self.disambiguated_grating = self.grating.lower()
        self.gratinglist = [self.grating]
        self.aperturelist = []
        self.instrumentlist = []
        self.propid = ''
        self.rootname = ''
        self.num_exp = 1

class Ullyses_COSSegmentList(COSSegmentList, Ullyses_SegmentList):
    pass


class Ullyses_STISSegmentList(STISSegmentList, Ullyses_SegmentList):
    pass


class Ullyses_CCDSegmentList(CCDSegmentList, Ullyses_SegmentList):
    pass


class Ullyses_FUSESegmentList(FUSESegmentList, Ullyses_SegmentList):
    pass


def find_files(indir): 
    allfiles = []
    for root, dirs, files in os.walk(indir, topdown=False):
        # Given a dir structure as follow, setting depth=2 ensure subdir/ will not be read
        # ULLYSES_DATA/
        # |___targ1/
        #     |___subdir/
        #
        depth = 2
        if root[len(indir):].count(os.sep) >= depth:
            continue

        print(root)
        dirname = root.split('/')[-1]
        print(f"   {dirname}")

        nonvofiles = glob.glob(os.path.join(root, '*_x1d.fits')) + glob.glob(os.path.join(root, '*_sx1.fits'))
        vofiles = glob.glob(os.path.join(root, '*_vo.fits'))
        allfiles += nonvofiles
        allfiles += vofiles
    return allfiles 


def coadd_and_abut_files(infiles, outdir, version=__release__, clobber=False):
    outdir_inplace = False
    if outdir is None:
        HLSP_DIR = os.getenv('HLSP_DIR')
        if HLSP_DIR is None:
            print("Environment variable HLSP_DIR must be defined if outdir is not specified")
            raise RuntimeError("Please set HLSP_DIR and restart")
        outdir_inplace = True
    
    nonvofiles = [x for x in infiles if "_vo.fits" not in x]
    vofiles = [x for x in infiles if "_vo.fits" in x]
    # collect the gratings that we will loop through
    # coadd.py will find the correct files itself,
    # but we need to know which gratings are present
    uniqmodes = defaultdict(list)
    
    for infile in nonvofiles:
        prihdr = fits.getheader(infile)
        obsmode = (prihdr['INSTRUME'], prihdr['OPT_ELEM'], prihdr['DETECTOR'])
        uniqmodes[obsmode].append(infile)

    if vofiles:
        if len(vofiles) != 1:
            print("More than 1 FUSE data file, aborting")
        else:
            obsmode = ('FUSE', 'FUSE', 'FUSE')
            uniqmodes[obsmode].append(vofiles[0])

    if not uniqmodes:
        print(f'No data to coadd, EXITING')
        return

    # Create dictionary of all products, with each set to None by default
    products = defaultdict(lambda: None)
    productdict = {}

    level = 2
    for obskey in uniqmodes:
        instrument, grating, detector = obskey
        infiles = uniqmodes[obskey]
        # this instantiates the class
        if instrument == 'COS':
            prod = Ullyses_COSSegmentList(grating, infiles=infiles)
        elif instrument == 'STIS':
            if detector == 'CCD':
                prod = Ullyses_CCDSegmentList(grating, infiles=infiles)
            else:
                prod = Ullyses_STISSegmentList(grating, infiles=infiles)
        elif instrument == 'FUSE':
            prod = Ullyses_FUSESegmentList(grating, infiles=infiles)
            products[f'{instrument}/{grating}'] = prod
        else:
            print(f'Unknown mode [{instrument}, {grating}, {detector}]')
            return


        # these two calls perform the main functions
        if len(prod.members) > 0:
            prod.create_output_wavelength_grid()
            prod.coadd()
            # this writes the output file
            # If making HLSPs for a DR, put them in the official folder
            prod.target = prod.get_targname()
            target = prod.target.lower()
            if outdir_inplace is True:
                outdir = os.path.join(HLSP_DIR, dir_target, version)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            if instrument != 'FUSE': # FUSE data is written as level 3 product below
                outname = create_output_file_name(prod, version, level=level)
                outname = os.path.join(outdir, outname)
                prod.write(outname, clobber, level=level, version=version)
                print(f"   Wrote {outname}")
            products[f'{instrument}/{grating}'] = prod
        else:
            print(f"No valid data for grating {grating}")
        if prod.level0 is True:
            prod.create_output_wavelength_grid()
            prod.coadd()
            # this writes the output file
            # If making HLSPs for a DR, put them in the official folder
            target = prod.target.lower()
            if outdir_inplace is True:
                outdir = os.path.join(HLSP_DIR, dir_target, version)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            outname = create_output_file_name(prod, version, level=0)
            outname = os.path.join(outdir, outname)
            prod.write(outname, clobber, level=0, version=version)
            print(f"   Wrote {outname}")
            products[f'{instrument}/{grating}'] = prod
        prod.add_hasp_attributes()
        products[f'{instrument}/{grating}'] = prod
        productdict[f'{instrument}/{grating}'] = prod

    # Create level 3 products- abutted spectra for gratings of the same
    # resolution for each instrument. Only create this product if more than
    # one grating of each resolution group is detected.
    level = 3
    lvl3_modes = {"cos_fuv_m": ["COS/G130M", "COS/G160M", "COS/G185M", "COS/G285M", "COS/G225M"],
                  "stis_m": ["STIS/E140M", "STIS/E230M", "STIS/G230MB", "STIS/G140M", "STIS/G230M", "STIS/G750M", "STIS/G430M"],
                  "stis_h": ["STIS/E140H", "STIS/E230H"],
                  "stis_l": ["STIS/G140L", "STIS/G230L", "STIS/G230LB", "STIS/G430L", "STIS/G750L"],
                  }
    for outprod,modes in lvl3_modes.items():
        lvl3_productlist = []
        lvl3_productdict = {}
        dowrite = False
        for mode in modes:
            if products[mode] is not None:
                lvl3_productlist.append(products[mode])
                lvl3_productdict[mode] = products[mode]
                dowrite = True
        # We only write level 3 products if more than one grating was found per resolution
        # group
        if dowrite is True and len(lvl3_productlist) > 1:
            lvl3_product = create_level4_products(lvl3_productlist, lvl3_productdict,
                                             grating_table=ULLYSES_GRATING_PRIORITIES)
            filename = create_output_file_name(lvl3_product, version, level=level)
            filename = os.path.join(outdir, filename)
            lvl3_product.write(filename, clobber, level=level, version=version)
            print(f"   Wrote {filename}")


    # Manually write out a FUSE level3 product.
    if products['FUSE/FUSE'] is not None:
        filename = create_output_file_name(products['FUSE/FUSE'], version, level=level)
        filename = os.path.join(outdir, filename)
        products['FUSE/FUSE'].write(filename, clobber, level=level, version=version)

    # Determine which gratings should contribute to the final level 4 SED HLSP.
    # Starting with the bluest product and working redward, find which products,
    # if any, overlap with the bluer product. If more than one overlaps, use
    # the one that extends further. If none overlap, still abut them- there
    # will just be a region of flux=0 in between.
    level = 4
    productlist = [productdict[key] for key in productdict]
    abutted_product = create_level4_products(productlist, productdict,
                                             grating_table=ULLYSES_GRATING_PRIORITIES)
    filename = create_output_file_name(abutted_product, version, level=level)
    filename = os.path.join(outdir, filename)
    abutted_product.write(filename, clobber, level=level, version=version)
    print(f"   Wrote {filename}")


def create_output_file_name(prod, version=__release__, level=3):
    instrument = prod.instrument.lower()   # will be either cos, stis, or fuse. If abbuted can be cos-stis or cos-stis-fuse
    grating = prod.grating.lower()
    target = prod.target.lower()
    version = version.lower()
    aperture = prod.aperture.lower()

    if level == 0:
        tel = 'hst'
        suffix = "spec"
    if level == 1:
        suffix = "mspec"
        tel = 'hst'
    elif level == 2:
        tel= 'hst'
        suffix = "cspec"
    elif level == 3:
        if instrument == 'fuse':
            tel = 'fuse'
            instrument = 'fuv'   # "fuv" is the "instrument" equivalent for fuse
            grating = aperture   # the grating for fuse data is set to "fuse" to change to use aperture
        else:
            tel= 'hst'
        suffix = 'aspec'
    elif level == 4:
        suffix = "preview-spec"
        if "g430" in prod.grating or "g750" in prod.grating:
            grating = "uv-opt"
        else:
            grating = "uv"
        if 'fuse' in instrument:
            tel = 'hst-fuse'
            # Move fuse to the front of the instrument list
            instrument_list = instrument.split('-')
            instrument_list.remove('fuse')
            instrument_list.insert(0, 'fuse')
            instrument = '-'.join(instrument_list)
        else:
            tel = 'hst'

    # Need to add logic for uv-opt here
    name = f"hlsp_ullyses_{tel}_{instrument}_{target}_{grating}_{version}_{suffix}.fits"
    return name


def main(indir, outdir, version=__release__, clobber=False):
    allfiles = find_files(indir)
    coadd_and_abut_files(allfiles, outdir, version, clobber)


def coadd_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",
                        default="./",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Directory for output HLSPs")
    parser.add_argument("-v", "--version", default=__release__,
                        help="Version number of the HLSP")
    parser.add_argument("-c", "--clobber", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")
    args = parser.parse_args()

    main(args.indir, args.outdir, args.version, args.clobber)


if __name__ == '__main__':
    coadd_parser()
