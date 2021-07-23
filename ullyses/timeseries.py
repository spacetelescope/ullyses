# Get a list of corrtag exposures with required grating
# Run splittag on each of these corrtags -> many new corrtags
# Run calcos on all of these corrtags
# - throw away everything other than the x1d files
# figure out the wavelength scale from all the x1ds
# coadd each of the x1d files into a new 1-d array of flux.  Each file will have the same wavelength array
# put the flux arrays into an array, one spectrum -> 1 row of array
import argparse
import astropy
from astropy.io import fits
from astropy.time import Time

import datetime
from datetime import datetime as dt
import re
import numpy as np
from . import coadd, wrapper
from costools import splittag
import calcos
import glob
import os

SECONDS_PER_DAY = 86400.0
version = wrapper.default_version

def collect_corrtag_inputs(grating):
    corrtag_list = []
    all_corrtags = glob.glob('split*_corrtag.fits')
    for file in all_corrtags:
        f1 = fits.open(file)
        if f1[0].header['OPT_ELEM'] == grating:
            corrtag_list.append(file)
        f1.close()
    return corrtag_list

def run_splittag(corrtag_list, duration):
    for file in corrtag_list:
        f1 = fits.open(file)
        rootname = f1[0].header['rootname']
        segment = f1[0].header['segment'][-1].lower()
        root = 'split_' + rootname
        splittag.splittag(file, root, starttime = 0., increment = duration, endtime = 1000.,
                          time_list = "")
    return

def run_calcos():
    splits = []
    all_splits = glob.glob('split*corrtag_a.fits')
    for split in all_splits:
        calcos.calcos(split, outdir='./x1dfiles')
        fltfiles = glob.glob('./x1dfiles/split*flt*.fits')
        for file in fltfiles:
            os.remove(file)
        countsfiles = glob.glob('./x1dfiles/split*counts*.fits')
        for file in countsfiles:
            os.remove(file)
    return

def sort_x1ds(min_exptime=20.0):
    x1dlist = []
    x1dfiles = glob.glob('split*_without.fits')
    for file in x1dfiles:
        f1 = fits.open(file)
        expend = f1[1].header['expend']
        exptime = f1[1].header['exptime']
        expstart = expend - exptime/SECONDS_PER_DAY
        if exptime > min_exptime:
            x1dlist.append((file, expstart))
        f1.close()
    sorted_x1dlist = sorted(x1dlist, key=lambda x: x[1])
    return sorted_x1dlist

def create_ensemble_segmentlist(grating, wavelength_binning=1):
    a = coadd.COSSegmentList(grating)
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
    return a

def rename_all_x1ds():
    x1dfiles = glob.glob('split*_x1d.fits')
    for oldfile in x1dfiles:
        newfile = oldfile.replace('_x1d.fits', '_without.fits')
        os.rename(oldfile, newfile)
    return

def transfer_from_ensemble(ensemble, segmentlist):
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

def process_all_files(grating, outfile, wavelength_binning=1):
    ensemble = create_ensemble_segmentlist(grating, wavelength_binning)
    rename_all_x1ds()
    withoutfiles = glob.glob('split*_without.fits')
    sorted_withoutfiles = sort_x1ds()
    nrows = len(sorted_withoutfiles)
    ncols = len(ensemble.output_wavelength)
    output_flux = np.zeros((nrows, ncols))
    output_error = np.zeros((nrows, ncols))
    row = 0
    starttimes = []
    endtimes = []
    wavelengths = ensemble.output_wavelength
    ensemble.primary_headers = []
    ensemble.first_headers = []
    for oldfile, expstart in sorted_withoutfiles:
        newfile = oldfile.replace('_without.fits', '_x1d.fits')
        os.rename(oldfile, newfile)
        a = coadd.COSSegmentList('G160M')
        transfer_from_ensemble(ensemble, a)
        a.coadd()
        start_time = a.first_headers[0]['expend'] - a.first_headers[0]['exptime']/SECONDS_PER_DAY
        endtimes.append(a.first_headers[0]['expend'])
        starttimes.append(start_time)
        output_flux[row] = a.output_flux
        output_error[row] = a.output_errors
        print(f'finished row {row}')
        row = row + 1
        os.rename(newfile, oldfile)
        ensemble.primary_headers.append(a.primary_headers[0])
        ensemble.first_headers.append(a.first_headers[0])
    for oldfile, expstart in sorted_withoutfiles:
        newfile = oldfile.replace('_without.fits', '_x1d.fits')
        os.rename(oldfile, newfile)
    write_product(output_flux, output_error, starttimes, endtimes, wavelengths, ensemble, outfile)
    return

def write_product(flux, error, starttimes, endtimes, wavelengths, ensemble, outfile):
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
    hdulist.writeto(outfile)
    return

def create_primary_header(ensemble, filename):
    level = 5
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
    hdr0['TARGNAME'] = ensemble.targname[0]
    hdr0.add_blank(after='OBSMODE')
    hdr0.add_blank('              / TARGET INFORMATION', before='TARGNAME')

    hdr0['RADESYS'] = ('ICRS ','World coordinate reference frame')
    ra, dec = ensemble.ull_coords()
    hdr0['TARG_RA'] =  (ra,  '[deg] Target right ascension')
    hdr0['TARG_DEC'] =  (dec,  '[deg] Target declination')
    hdr0['PROPOSID'] = (ensemble.combine_keys("proposid", "multi"), 'Program identifier')
    hdr0.add_blank(after='TARG_DEC')
    hdr0.add_blank('           / PROVENANCE INFORMATION', before='PROPOSID')
    hdr0['CAL_VER'] = (f'ULLYSES Cal {coadd.CAL_VER}', 'HLSP processing software version')
    hdr0['HLSPID'] = ('ULLYSES', 'Name ID of this HLSP collection')
    hdr0['HSLPNAME'] = ('Hubble UV Legacy Library of Young Stars as Essential Standards',
                    'Name ID of this HLSP collection')
    hdr0['HLSPLEAD'] = ('Julia Roman-Duval', 'Full name of HLSP project lead') 
    hdr0['HLSP_VER'] = (version,'HLSP data release version identifier')
    hdr0['HLSP_LVL'] = (level, 'ULLYSES HLSP Level')
    hdr0['LICENSE'] = ('CC BY 4.0', 'License for use of these data')
    hdr0['LICENURL'] = ('https://creativecommons.org/licenses/by/4.0/', 'Data license URL')
    hdr0['REFERENC'] = ('https://ui.adsabs.harvard.edu/abs/2020RNAAS...4..205R', 'Bibliographic ID of primary paper')
    
    hdr0['CENTRWV'] = (ensemble.combine_keys("centrwv", "average"), 'Central wavelength of the data')
    hdr0.add_blank(after='REFERENC')
    hdr0.add_blank('           / ARCHIVE SEARCH KEYWORDS', before='CENTRWV')
    hdr0['MINWAVE'] = (ensemble.combine_keys("minwave", "min"), 'Minimum wavelength in spectrum')
    hdr0['MAXWAVE'] = (ensemble.combine_keys("maxwave", "max"), 'Maximum wavelength in spectrum')

    primary = fits.PrimaryHDU(header=hdr0)
    return primary
    
def create_extension_1_header(ensemble):
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
    return hdr1
    
def create_extension_2(ensemble):
    hdr2 = fits.Header()
    hdr2['EXTNAME'] = ('PROVENANCE', 'Metadata for contributing observations')
    # set up the table columns
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
    mjd_begs = ensemble.combine_keys("expstart", "arr")
    mjd_ends = ensemble.combine_keys("expend", "arr")
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

def main(run_calcos, wavelength_binning, grating, resolution, normal, outfile):
    
    if run_calcos:
        corrtag_list = collect_corrtag_inputs(grating)
        duration = resolution
        run_splittag(corrtag_list, duration)
        run_calcos()
    process_all_files(grating, outfile, wavelength_binning=wavelength_binning)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--run_calcos", default=True,
                        help="Run splittag and calcos?")
    parser.add_argument("-b", "--wavelength_binning", default=1.,
                        help="Wavelength binning")
    parser.add_argument("-g", "--grating", default='',
                        help="Grating")
    parser.add_argument("-r", "--resolution", default=30.0, 
                        help="Time resolution in seconds for timeseries product")
    parser.add_argument("-n", "--normal", default=True,
                        help="Create normal coadd product from all inputs?")
    parser.add_argument("-o", "--output", default='',
                        help="Name of output file for timeseries product")
    args = parser.parse_args()

    main(args.run_calcos, args.wavelength_binning, args.grating, args.resolution, args.normal,
         args.output)
