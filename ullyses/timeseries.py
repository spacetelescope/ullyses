# Get a list of corrtag exposures with required grating
# Run splittag on each of these corrtags -> many new corrtags
# Run calcos on all of these corrtags
# - throw away averything other than the x1d files
# figure out the wavelength scale from all the x1ds
# coadd each of the x1d files into a new 1-d array of flux.  Each file will have the same wavelength array
# put the flux arrays into an array, one spectrum -> 1 row of array

from astropy.io import fits
import numpy as np
from high_level_science_products import coadd
from costools import splittag
import calcos
import glob
import os

SECONDS_PER_DAY = 86400.0

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

def create_ensemble_segmentlist(wavelength_binning=1):
    a = coadd.COSSegmentList('G160M')
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

def process_all_files(outfile, wavelength_binning=1):
    ensemble = create_ensemble_segmentlist(wavelength_binning)
    rename_all_x1ds()
    withoutfiles = glob.glob('split*_without.fits')
    sorted_withoutfiles = sort_x1ds()
    nrows = len(sorted_withoutfiles)
    ncols = len(ensemble.output_wavelength)
    output_flux = np.zeros((nrows, ncols))
    output_error = np.zeros((nrows, ncols))
    row = 0
    times = []
    wavelengths = ensemble.output_wavelength
    for oldfile, expstart in sorted_withoutfiles:
        newfile = oldfile.replace('_without.fits', '_x1d.fits')
        os.rename(oldfile, newfile)
        a = coadd.COSSegmentList('G160M')
        transfer_from_ensemble(ensemble, a)
        a.coadd()
        start_time = a.first_headers[0]['expend'] - a.first_headers[0]['exptime']/SECONDS_PER_DAY
        times.append(start_time)
        output_flux[row] = a.output_flux
        output_error[row] = a.output_errors
        print(f'finished row {row}')
        row = row + 1
        os.rename(newfile, oldfile)
    for oldfile, expstart in sorted_withoutfiles:
        newfile = oldfile.replace('_without.fits', '_x1d.fits')
        os.rename(oldfile, newfile)
    write_product(output_flux, output_error, times, wavelengths, outfile)
    return

def write_product(flux, error, times, wavelengths, outfile):
    nrows, ncolumns = flux.shape
    npixels = nrows*ncolumns
    columns = []
    array_dimensions = '(' + str(ncolumns) + ', ' + str(nrows) + ')'
    columns.append(fits.Column(name='TIME', format=str(nrows)+'D', unit='Day'))
    columns.append(fits.Column(name='WAVELENGTH', format=str(ncolumns)+'E', unit='Angstrom'))
    columns.append(fits.Column(name='FLUX', format=str(npixels)+'E', dim=array_dimensions, unit='erg /s /cm**2 /Angstrom'))
    columns.append(fits.Column(name='ERROR', format=str(npixels)+'E', dim=array_dimensions, unit='erg /s /cm**2 /Angstrom'))
    cd = fits.ColDefs(columns)
    table1 = fits.BinTableHDU.from_columns(cd, nrows=1)
    table1.data[0]['wavelength'] = wavelengths
    table1.data[0]['time'] = times
    table1.data[0]['flux'] = flux
    table1.data[0]['error'] = error
    prihdu = fits.PrimaryHDU()
    hdulist = fits.HDUList()
    hdulist.append(prihdu)
    hdulist.append(table1)
    hdulist.writeto(outfile)
    return

if __name__ == '__main__':
    grating = 'G160M'
    corrtag_list = collect_corrtag_inputs(grating)
    duration = 50.0
    run_splittag(corrtag_list, duration)
    run_calcos()

    
