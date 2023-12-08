import argparse
import astropy
from astropy.io import fits
from astropy.time import Time
import datetime
from datetime import datetime as dt
import re
import numpy as np
from costools import splittag
import calcos
import glob
import os

from ullyses import coadd
from ullyses import ullyses_coadd_abut_wrapper as wrapper
from . import __version__, __release__

SECONDS_PER_DAY = 86400.0
RED = "\033[1;31m"
RESET = "\033[0;0m"

"""
    Timeseries product creation code.

    This code should be run from a directory that contains only files from 1
    grating, and either split or full files (but not both).  A suitable set of
    directories to run this from are:

    base_directory/G160M/full/
    base_directory/G160M/split/
    base_directory/G230L/full/
    base_directory/G230L/split/

    Put the x1d files that you want to contribute to the product into the appropriate
    directories.  If you want to exclude files, leave them out from the directory.

    You will need to install this code using python setup.py install

    To create a timeseries product, go to the selected directory

    Start python/ipython

>>> from high_level_science_products import timeseries
>>> timeseries.process_files('G160M', 'timeseries_product.fits', wavelength_binning=3, min_exptime=20.)

    That's all that should be necessary to create the timeseries product named
    'timeseries_product.fits' from the x1d files in the current directory with the
    specified grating and with exposure times greater than specified in the call

"""

def sort_split_x1ds(ins, grating, indir=".", min_exptime=20.0):
    """
    THIS IS NOT CURRENTLY IN USE.
    Sort the split x1d files.    Select all files that match the
    pattern 'split_*_without.fits' and that match the grating and whose
    exposure time is greated than min_exptime, and sort them by expstart

    Parameters:
    -----------

    grating: str
        Grating that must be present in OPT_ELEM keyword of primary
        header of each file

    min_exptime: float, default=20.0
        Minimum exposure time to be included

    Returns:
        list of tuples (filename, expstart), sorted by expstart

    """
    x1dfiles = glob.glob(os.path.join(indir, 'split*_without.fits'))
    good_list = []
    for file in x1dfiles:
        this_grating = fits.getval(file, 'OPT_ELEM')
        this_ins = fits.getval(file, "INSTRUME")
        if this_grating == grating and this_ins == ins:
            good_list.append(file)
    sorted_x1dlist = sort_x1dfiles(good_list, min_exptime=min_exptime)
    return sorted_x1dlist

def sort_full_x1ds(ins, grating, indir=".", min_exptime=20.0):
    """
    THIS IS NOT CURRENTLY IN USE.
    Sort the full x1d files.    Select all files that match the
    pattern '*_without.fits' and don't start with 'split', and that match
    the grating and whose exposure time is greated than min_exptime, and
    sort them by expstart

    Parameters:
    -----------

    grating: str
        Grating that must be present in OPT_ELEM keyword of primary
        header of each file

    min_exptime: float, default=20.0
        Minimum exposure time to be included

    Returns:
        list of tuples (filename, expstart), sorted by expstart

    """
    x1dfiles = glob.glob(os.path.join(indir, '*_without.fits'))
    good_list = []
    for file in x1dfiles:
        if file.startswith('split'):
            continue
        this_grating = fits.getval(file, 'OPT_ELEM')
        this_ins = fits.getval(file, "INSTRUME")
        if this_grating == grating and this_ins == ins:
            good_list.append(file)
    sorted_x1dlist = sort_x1dfiles(good_list, min_exptime=min_exptime)
    return sorted_x1dlist

def sort_x1ds(ins, grating, indir=".", min_exptime=20.0):
    """Sort the x1d files.    Select all files that match the
    pattern '*_without.fits' and that match
    the grating and whose exposure time is greater than min_exptime, and
    sort them by expstart

    Parameters:
    -----------

    grating: str
        Grating that must be present in OPT_ELEM keyword of primary
        header of each file

    min_exptime: float, default=20.0
        Minimum exposure time to be included

    Returns:
        list of tuples (filename, expstart), sorted by expstart

    """
    x1dfiles = glob.glob(os.path.join(indir, '*_without.fits'))
    good_list = []
    for file in x1dfiles:
        this_grating = fits.getval(file, 'OPT_ELEM')
        this_ins = fits.getval(file, "INSTRUME")
        if this_grating == grating and this_ins == ins:
            good_list.append(file)
    sorted_x1dlist = sort_x1dfiles(good_list, min_exptime=min_exptime)
    return sorted_x1dlist


def sort_x1dfiles(x1dfiles, min_exptime=20.0):
    """Sort a list of files in order of expstart

    Parameters:
    -----------

    x1dfiles: list
        List of filenames to be sorted

    min_exptime: float, default=20.0
        Minimum exposure time to be included in sorted output

    Returns:
        List of (filename, exptime) tuples

    """
    x1dlist = []
    for file in x1dfiles:
        expend = fits.getval(file, 'expend', 1)
        exptime = fits.getval(file, 'exptime', 1)
        expstart = expend - exptime/SECONDS_PER_DAY
        if exptime > min_exptime:
            x1dlist.append((file, expstart))
    sorted_x1dlist = sorted(x1dlist, key=lambda x: x[1])
    return sorted_x1dlist

def create_ensemble_segmentlist(grating, indir=".", wavelength_binning=1.0, ins="COS", infiles=None):
    """Create the ensemble Ullyses_COSSegmentList whose wavelength and array
    sizing parameters are to be used to create the output product.
    Should be created from full exposures (not splits), so make sure
    all split exposures are renamed from ending with _x1d.fits so they
    don't add to the Ullyses_COSSegmentList

    Parameters:
    -----------

    grating: str
        Grating to be used to create the timeseries product

    wavelength_binning: float, default=1.0
        Wavelength binning for output product in pixels.  Can be non-integer

    Returns:
        Ullyses_COSSegmentList:

    """
    if ins == "COS":
        a = wrapper.Ullyses_COSSegmentList(grating, inpath=indir, infiles=infiles)
    elif ins == "STIS":
        if grating in ["G230LB", "G230M", "G230LB", "G230MB", "G430L", "G430M", "G750L", "G750M"]:
            a = wrapper.Ullyses_CCDSegmentList(grating, inpath=indir, infiles=infiles)
        else:
            a = wrapper.Ullyses_STISSegmentList(grating, inpath=indir, infiles=infiles)
    a.create_output_wavelength_grid()
    if wavelength_binning != 1:
        a.delta_wavelength = a.delta_wavelength * wavelength_binning
        a.output_wavelength = np.arange(a.min_wavelength, a.max_wavelength, a.delta_wavelength)
        nelements = len(a.output_wavelength)
        a.output_sumgcounts = np.zeros(nelements)
        a.output_sumflux = np.zeros(nelements)
        a.output_sumweight = np.zeros(nelements)
        a.output_varsum = np.zeros(nelements)
        a.output_flux = np.zeros(nelements)
        a.output_errors = np.zeros(nelements)
        a.signal_to_noise = np.zeros(nelements)
        a.sumnetcounts = np.zeros(nelements)
        a.output_exptime = np.zeros(nelements)
    # sort the headers in order of expstart

    return a

def rename_all_split_x1ds(indir="."):
    """Rename all split exposures (that start with 'split') from
    ending with '_x1d.fits' to ending with '_without.fits'.  This stops
    the files from getting added to Ullyses_COSSegmentLists

    Parameters:
    -----------

    None

    Returns:
    --------

    None

    """
    x1dfiles = glob.glob(os.path.join(indir, 'split*_x1d.fits')) + glob.glob(os.path.join(indir, 'split*_sx1.fits'))
    for oldfile in x1dfiles:
        if "x1d.fits" in oldfile:
            newfile = oldfile.replace('_x1d.fits', '_without.fits')
        elif "sx1.fits" in oldfile:
            newfile = oldfile.replace('_sx1.fits', '_without.fits')
        os.rename(oldfile, newfile)
    return

def rename_all_x1ds(indir="."):
    """Rename all exposures from
    ending with '_x1d.fits' to ending with '_without.fits'.  This stops
    the files from getting added to Ullyses_COSSegmentLists

    Parameters:
    -----------

    None

    Returns:
    --------

    None

    """
    x1dfiles = glob.glob(os.path.join(indir, '*_x1d.fits')) + glob.glob(os.path.join(indir, '*_sx1.fits'))
    renaming_old = {}
    renaming_new = {}
    for oldfile in x1dfiles:
        if "x1d.fits" in oldfile:
            newfile = oldfile.replace('_x1d.fits', '_without.fits')
        elif "sx1.fits" in oldfile:
            newfile = oldfile.replace('_sx1.fits', '_without.fits')
        renaming_old[oldfile] = newfile
        renaming_new[newfile] = oldfile
        os.rename(oldfile, newfile)
    return renaming_old, renaming_new

def rename_all_full_x1ds(indir="."):
    """Rename all full exposures (that don't start with 'split') from
    ending with '_x1d.fits' to ending with '_without.fits'.  This stops
    the files from getting added to Ullyses_COSSegmentLists

    Parameters:
    -----------

    None

    Returns:
    --------

    None

    """
    x1dfiles = glob.glob(os.path.join(indir, '*_x1d.fits')) + glob.glob(os.path.join(indir, '*_sx1.fits'))
    for file in x1dfiles:
        if file.startswith('split'):
            x1dfiles.remove(file)
    for oldfile in x1dfiles:
        if "x1d.fits" in oldfile:
            newfile = oldfile.replace('_x1d.fits', '_without.fits')
        elif "sx1.fits" in oldfile:
            newfile = oldfile.replace('_sx1.fits', '_without.fits')
        os.rename(oldfile, newfile)
    return

def transfer_from_ensemble(ensemble, segmentlist):
    """Transfer the wavelength and array size information from the ensemble
    Ullyses_COSSegmentList to the target Ullyses_COSSegmentList

    Parameters:
    -----------

    ensemble: Ullyses_COSSegmentList
        Ullyses_COSSegmentList containing the wavelength and array sizing parameters to be used
        in the target Ullyses_COSSegmentList

    segmentlist: Ullyses_COSSegmentList
        Target Ullyses_COSSegmentList to be changed by transferring parameters from ensemble

    Returns:
    --------

    None

    """
    segmentlist.delta_wavelength = ensemble.delta_wavelength
    segmentlist.min_wavelength = ensemble.min_wavelength
    segmentlist.max_wavelength = ensemble.max_wavelength
    segmentlist.output_wavelength = ensemble.output_wavelength
    segmentlist.nelements = len(segmentlist.output_wavelength)
    segmentlist.output_sumgcounts = np.zeros(segmentlist.nelements)
    segmentlist.output_sumflux = np.zeros(segmentlist.nelements)
    segmentlist.output_sumweight = np.zeros(segmentlist.nelements)
    segmentlist.output_varsum = np.zeros(segmentlist.nelements)
    segmentlist.output_flux = np.zeros(segmentlist.nelements)
    segmentlist.output_errors = np.zeros(segmentlist.nelements)
    segmentlist.signal_to_noise = np.zeros(segmentlist.nelements)
    segmentlist.sumnetcounts = np.zeros(segmentlist.nelements)
    segmentlist.output_exptime = np.zeros(segmentlist.nelements)
    return

def process_files(grating, outfile, indir=".", wavelength_binning=1, min_exptime=20,
                  overwrite=False, ins="COS", infiles=None):
    """Process all x1d files in the current directory with the selected grating.

    Parameters:
    -----------

    grating: str
        The grating to be processed.  This should be the same as the OPT_ELEM keyword in
        the primary header of the files to be processed

    outfile: str
        The name of the output file for the timeseries product made from the x1d exposures

    wavelength_binning: float, default = 1.0
        The wavelength binning in pixels, can be a non-integer

    min_exptime: float, default=20.0
        The minimum exposure time for the splittag exposures.  Some of the splittag files
        from the end of an exposure don't have the full exposure time and are very noisy.
        This parameter excludes all files with exposure times that are below the specified
        threshold

    """
    # Create a Ullyses_COSSegmentList from all files that end with _x1d.fits and
    # have the appropriate grating. The ensemble
    # is used to determine the start, stop and delta wavelength of the product
    # and to contain the headers used to create the provenance table
    # Ullyses_COSSegmentLists have the Ullyses-specific get_target() and get_coords()
    # methods
    ensemble = create_ensemble_segmentlist(grating, indir, wavelength_binning, ins=ins, infiles=infiles)
    # If there's only 1 matching file for hte specified grating, exit 
    if len(ensemble.datasets) == 1:
        print(f"{RED}WARNING: only one matching dataset for grating {grating}, skipping{RESET}")
        return
    # Rename all the x1d.fits files to remove the _x1d.fits ending so that
    # subsequent Ullyses_COSSegmentLists can be created one at a time
    renaming_old, renaming_new = rename_all_x1ds(indir)
    # create a list of split files sorted by expstart
    sorted_files = sort_x1ds(ins, grating, indir, min_exptime=min_exptime)
    # process this list of files one at a time to each write a single row
    # to the output flux and error arrays.  Write the output file when complete
    process_sorted_filelist(sorted_files, grating, ensemble, outfile, renaming_new, indir,
                            overwrite, ins=ins)
    # Rename the files back from _without.fits to _x1d.fits
    rename_files_back(renaming_new)

def process_sorted_filelist(sorted_list, grating, ensemble, outfile, renaming_new, indir=".",
                            overwrite=False, ins="COS"):
    """Create a timeseries product from a sorted list of (input file, expstart) tuples.

    Parameters:
    -----------

    sorted_list: list of tuples
        List of (filename, expstart) tuples sorted by expstart.

    grating: str
        Grating to be processed.  All the input files SHOULD have the same grating,
        but just in case this allows rejection of exposures that don't match

    ensemble: Ullyses_COSSegmentList
        Ullyses_COSSegmentList made from exposures that are used to define the output wavelength
        grid.  The start wavelength, stop wavelength and wavelength increment are used to
        initialize and size output arrays

    outfile: str
        Name of output FITS file to write the product to

    Returns:

    None

    """
    nrows = len(sorted_list)
    ncols = len(ensemble.output_wavelength)
    output_flux = np.zeros((nrows, ncols))
    output_error = np.zeros((nrows, ncols))
    row = 0
    starttimes = []
    endtimes = []
    wavelengths = ensemble.output_wavelength
    ensemble.primary_headers = []
    ensemble.first_headers = []
    if ins == "COS":
        segmentlist = wrapper.Ullyses_COSSegmentList
    elif ins == "STIS":
        if grating in ["G230LB", "G230M", "G230LB", "G230MB", "G430L", "G430M", "G750L", "G750M"]:
            segmentlist = wrapper.Ullyses_CCDSegmentList
        else:
            segmentlist = wrapper.Ullyses_STISSegmentList
    for newfile, expstart in sorted_list:
        oldfile = renaming_new[newfile]
        os.rename(newfile, oldfile)
        # Make sure this file uses the required grating
        this_grating = fits.getval(oldfile, 'OPT_ELEM')
        this_ins = fits.getval(oldfile, "INSTRUME")
        if this_grating != grating and this_ins != ins:
            print(f"Skipping file {newfile} as it doesn't have the required grating")
            continue
        a = segmentlist(grating, inpath=indir)
        transfer_from_ensemble(ensemble, a)
        a.coadd()
        start_time = a.first_headers[0]['expend'] - a.first_headers[0]['exptime']/SECONDS_PER_DAY
        endtimes.append(a.first_headers[0]['expend'])
        starttimes.append(start_time)
        output_flux[row] = a.output_flux
        output_error[row] = a.output_errors
        print(f'finished row {row}')
        row = row + 1
        # Rename the file back to _without.fits so that it won't appear in subsequent
        # Ullyses_COSSegmentLists
        os.rename(oldfile, newfile)
        ensemble.primary_headers.append(a.primary_headers[0])
        ensemble.first_headers.append(a.first_headers[0])
    write_product(output_flux, output_error, starttimes, endtimes, wavelengths,
                  ensemble, outfile, overwrite)
    return

def rename_files_back(renaming_new):
    """After processing files, they are left with names that end in '_without.fits'.
    This routine renames them back so they end in '_x1d.fits'

    Parameters:
    -----------

    None

    Returns:
    --------

    None

    """
    for newfile,oldfile in renaming_new.items():
        os.rename(newfile, oldfile)

def write_product(flux, error, starttimes, endtimes, wavelengths, ensemble, outfile, overwrite=False):
    """Write the timeseries product to a FITS file.  The output file has a primary extension
    with no data, and 2 data extensions.  The first data extension contains a table with 5
    columns: EXPSTART, EXPEND, WAVELENGTH, FLUX, ERROR.  EXPSTART and EXPEND give the start
    and end times of each row of the FLUX and ERROR arrays, while WAVELENGTH gives the
    wavelength of each column of these arrays.  The second extension contains the provenance
    information about the files that were used to create these products.

    Parameters:
    -----------

    flux: float ndarray
        2-d array of flux values

    error: float ndarray
        2-d array of error values

    starttimes: float ndarray
        1-d array of start times (MJD)

    endtimes: float ndarray
        1-d array of stop times (MJD)

    wavelengths: float ndarray
        1-d array of wavelengths

    ensemble: Ullyses_COSSegmentList
        Ullyses_COSSegmentList made from exposures that were used to create this product

    outfile: str
        Name of FITS product file

    Returns:
    --------

    None

    """
    nrows, ncolumns = flux.shape
    npixels = nrows*ncolumns
    columns = []
    array_dimensions = '(' + str(ncolumns) + ', ' + str(nrows) + ')'
    columns.append(fits.Column(name='MJDSTART', format=str(nrows)+'D', unit='Day'))
    columns.append(fits.Column(name='MJDEND', format=str(nrows) + 'D', unit='Day'))
    columns.append(fits.Column(name='WAVELENGTH', format=str(ncolumns)+'E', unit='Angstrom'))
    columns.append(fits.Column(name='FLUX', format=str(npixels)+'E', dim=array_dimensions, unit='erg /s /cm**2 /Angstrom'))
    columns.append(fits.Column(name='ERROR', format=str(npixels)+'E', dim=array_dimensions, unit='erg /s /cm**2 /Angstrom'))
    cd = fits.ColDefs(columns)
    hdr1 = create_extension_1_header(ensemble)
    table1 = fits.BinTableHDU.from_columns(cd, header=hdr1, nrows=1)
    table1.data[0]['wavelength'] = wavelengths
    table1.data[0]['mjdstart'] = starttimes
    table1.data[0]['mjdend'] = endtimes
    table1.data[0]['flux'] = flux
    table1.data[0]['error'] = error
    prihdu = create_primary_header(ensemble, outfile)
    hdulist = fits.HDUList()
    hdulist.append(prihdu)
    hdulist.append(table1)
    table2 = create_extension_2(ensemble)
    hdulist.append(table2)
    hdulist.writeto(outfile, overwrite=overwrite)

    print(f"\nWrote output file: {outfile}")
    return

def create_primary_header(ensemble, filename):
    """Create the primary header of the timeseries product.

    Parameters:

    ensemble: Ullyses_COSSegmentList
        Ullyses_COSSegmentList of full exposures used to make the arrays in the timeseries
        product.  The primary_headers and first_headers attributes are used to
        populate the header

    filename: str
        Name of timeseries product filename

    Returns:

    prihdr: astropy.io,fits PrimaryHDU object

    """
    level = 5

    # If the target is a ULLYSES target, use the official
    # target name and coords
    ensemble.target = ensemble.get_targname()
    ensemble.targ_ra, ensemble.targ_dec, ensemble.coord_epoch = ensemble.get_coords()

    hdr0 = fits.Header()
    hdr0['EXTEND'] = ('T', 'FITS file may contain extensions')
    hdr0['NEXTEND'] = 3
    hdr0['FITS_VER'] = 'Definition of the Flexible Image Transport System (FITS) v4.0 https://fits.gsfc.nasa.gov/standard40/fits_standard40aa-le.pdf'
    hdr0['FITS_SW'] = ('astropy.io.fits v' + astropy.__version__, 'FITS file creation software')
    hdr0['ORIGIN'] = ('Space Telescope Science Institute', 'FITS file originator')
    hdr0['DATE'] = (str(datetime.date.today()), 'Date this file was written')
    hdr0['FILENAME'] = (os.path.basename(filename), 'Name of this file')
    hdr0['TELESCOP'] = (ensemble.combine_keys("telescop", "multi"), 'Telescope used to acquire data')
    hdr0['INSTRUME'] = (ensemble.combine_keys("instrume", "multi"), 'Instrument used to acquire data')
    hdr0.add_blank('', after='TELESCOP')
    hdr0.add_blank('              / SCIENCE INSTRUMENT CONFIGURATION', before='INSTRUME')
    hdr0['DETECTOR'] = (ensemble.combine_keys("detector", "multi"), 'Detector or channel used to acquire data')
    hdr0['DISPERSR'] = (ensemble.combine_keys("opt_elem", "multi"), 'Identifier of disperser')
    hdr0['CENWAVE'] = (ensemble.combine_keys("cenwave", "multi"), 'Central wavelength setting for disperser')
    hdr0['APERTURE'] = (ensemble.combine_keys("aperture", "multi"), 'Identifier of entrance aperture')
    hdr0['S_REGION'] = (ensemble.obs_footprint(), 'Region footprint')
    hdr0['OBSMODE'] = (ensemble.combine_keys("obsmode", "multi"), 'Instrument operating mode (ACCUM | TIME-TAG)')
    hdr0['TARGNAME'] = ensemble.target
    hdr0.add_blank(after='OBSMODE')
    hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

    hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
    ra, dec, coord_epoch = ensemble.get_coords()
    hdr0['G_EPOCH'] =  (coord_epoch,  'Epoch of GAIA coordinates')
    hdr0['TARG_RA'] =  (ra,  '[deg] Target right ascension')
    hdr0['TARG_DEC'] =  (dec,  '[deg] Target declination')
    hdr0['PROPOSID'] = (ensemble.combine_keys("proposid", "multi"), 'Program identifier')
    hdr0.add_blank(after='TARG_DEC')
    hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
    hdr0['CAL_VER'] = (f'ULLYSES Cal {__version__}', 'HLSP processing software version')
    hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
    hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                    'Name ID of this HLSP collection')
    hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead')
    hdr0['HLSP_VER'] = (__release__,'HLSP data release version identifier')
    hdr0['HLSP_LVL'] = (level, 'ULLYSES HLSP Level')
    hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
    hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
    hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')

    minwave = ensemble.combine_keys("minwave", "min")
    maxwave = ensemble.combine_keys("maxwave", "max")
    centrwv = ((maxwave - minwave)/2.) + minwave
    hdr0['CENTRWV'] = (centrwv, 'Central wavelength of the data')
    hdr0.add_blank(after='REFERENC')
    hdr0.add_blank('           / ARCHIVE SEARCH KEYWORDS', before='CENTRWV')
    hdr0['MINWAVE'] = (minwave, 'Minimum wavelength in spectrum')
    hdr0['MAXWAVE'] = (maxwave, 'Maximum wavelength in spectrum')

    primary = fits.PrimaryHDU(header=hdr0)
    return primary

def create_extension_1_header(ensemble):
    """Create the extension 1 header for the timeseries product

    Parameters:
    -----------

    ensemble: Ullyses_COSSegmentList
        The Ullyses_COSSegmentList from which the header parameters are derived.  This should be
        created from the full exposures

    Returns:
    --------

    astropy.io.fits.Header object

    """
    hdr1 = fits.Header()
    hdr1['EXTNAME'] = ('SCIENCE', 'Spectrum science arrays')
    hdr1['TIMESYS'] = ('UTC', 'Time system in use')
    hdr1['TIMEUNIT'] = ('s', 'Time unit for durations')
    hdr1['TREFPOS'] = ('GEOCENTER', 'Time reference position')

    mjd_beg = ensemble.combine_keys("expstart", "min")
    mjd_end = ensemble.combine_keys("expend", "max")
    dt_beg = Time(mjd_beg, format="mjd").datetime
    dt_end = Time(mjd_end, format="mjd").datetime
    hdr1['DATE-BEG'] = (dt.strftime(dt_beg, "%Y-%m-%dT%H:%M:%S"), 'Date-time of first observation start')
    hdr1.add_blank('', after='TREFPOS')
    hdr1.add_blank('              / FITS TIME COORDINATE KEYWORDS', before='DATE-BEG')

    hdr1['DATE-END'] = (dt.strftime(dt_end, "%Y-%m-%dT%H:%M:%S"), 'Date-time of last observation end')
    hdr1['MJD-BEG'] = (mjd_beg, 'MJD of first exposure start')
    hdr1['MJD-END'] = (mjd_end, 'MJD of last exposure end')
    hdr1['XPOSURE'] = (ensemble.combine_keys("exptime", "sum"), '[s] Sum of exposure durations')
    all_comments = ensemble.combine_keys("comment", "comment")
    for comment in all_comments:
        hdr1['COMMENT'] = (comment, "Calibration and/or quality comment")
    return hdr1

def create_extension_2(ensemble):
    """Create the timeseries product provenance extension

    Parameters:
    -----------

    ensemble: Ullyses_COSSegmentList
        The Ullyses_COSSegmentList from which the provenance data are to be derived.  Should
        be created from the full exposures.

    Returns:

    astropy.io.fits.BinTableHDU object

    """
    hdr2 = fits.Header()
    hdr2['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
    # set up the table columns
    # If this product is created from the splittag files, the FILENAME column needs to be updated
    # so that it points to the corrtag files that were used to create the timeseries product
    first_filename = ensemble.primary_headers[0]['FILENAME']
    if first_filename.startswith('split'):
        detector = ensemble.primary_headers[0]['DETECTOR']
        filename_array = ensemble.combine_keys("filename", "arr")
        new_filename_list = []
        for x1dfilename in filename_array:
            ipppssoot = x1dfilename[6:15]
            if detector == 'FUV':
                replacement_string = '_corrtag_a.fits'
            else:
                replacement_string = '_corrtag.fits'
            newfilename = ipppssoot + replacement_string
            new_filename_list.append(newfilename)
        filename_array = np.array(new_filename_list)
        cfn = fits.Column(name='FILENAME', array=filename_array, format='A40')
    else:
        cfn = fits.Column(name='FILENAME', array=ensemble.combine_keys("filename", "arr"), format='A40')
    cpid = fits.Column(name='PROPOSID', array=ensemble.combine_keys("proposid", "arr"), format='A32')
    ctel = fits.Column(name='TELESCOPE', array=ensemble.combine_keys("telescop", "arr"), format='A32')
    cins = fits.Column(name='INSTRUMENT', array=ensemble.combine_keys("instrume", "arr"), format='A32')
    cdet = fits.Column(name='DETECTOR', array=ensemble.combine_keys("detector", "arr"), format='A32')
    cdis = fits.Column(name='DISPERSER', array=ensemble.combine_keys("opt_elem", "arr"), format='A32')
    ccen = fits.Column(name='CENWAVE', array=ensemble.combine_keys("cenwave", "arr"), format='A32')
    cap = fits.Column(name='APERTURE', array=ensemble.combine_keys("aperture", "arr"), format='A32')
    csr = fits.Column(name='SPECRES', array=ensemble.combine_keys("specres", "arr"), format='F8.1')
    ccv = fits.Column(name='CAL_VER', array=ensemble.combine_keys("cal_ver", "arr"), format='A32')
    mjd_ends = ensemble.combine_keys("expend", "arr")
    exptime = ensemble.combine_keys("exptime", "arr")
    mjd_begs = mjd_ends - exptime / SECONDS_PER_DAY
    mjd_mids = (mjd_ends + mjd_begs) / 2.
    cdb = fits.Column(name='MJD_BEG', array=mjd_begs, format='F15.9', unit='d')
    cdm = fits.Column(name='MJD_MID', array=mjd_mids, format='F15.9', unit='d')
    cde = fits.Column(name='MJD_END', array=mjd_ends, format='F15.9', unit='d')
    cexp = fits.Column(name='XPOSURE', array=ensemble.combine_keys("exptime", "arr"), format='F15.9', unit='s')
    cmin = fits.Column(name='MINWAVE', array=ensemble.combine_keys("minwave", "arr"), format='F9.4', unit='Angstrom')
    cmax = fits.Column(name='MAXWAVE', array=ensemble.combine_keys("maxwave", "arr"), format='F9.4', unit='Angstrom')

    cd2 = fits.ColDefs([cfn, cpid, ctel, cins, cdet, cdis, ccen, cap, csr, ccv, cdb, cdm, cde, cexp, cmin ,cmax])

    table2 = fits.BinTableHDU.from_columns(cd2, header=hdr2)
    return table2


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default=".",
                        help="Path to input directory with either default x1d files or split x1d files")
    parser.add_argument("-g", "--grating", 
                        help="Grating to process")
    parser.add_argument("-o", "--outfile",
                        help="Name of output file")
    parser.add_argument("-w", "--wl", default=1, type=int,
                        help="Number of pixels to bin in X")
    parser.add_argument("-e", "--min_exp", default=20, type=int,
                        help="Minimum required exposure time of split corrtag")
    parser.add_argument("-c", "--clobber", default=False, action="store_true",
                        help="If True, overwrite existing products")
    args = parser.parse_args()

    process_files(args.grating, args.outfile, args.indir, args.wl, args.min_exp,
         args.clobber)
