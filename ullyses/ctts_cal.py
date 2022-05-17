import argparse
import numpy as np
from astropy.io import fits
import glob
import os
import shutil
import calcos

from ullyses import splittag_wrapper
from ullyses import timeseries
import ullyses_utils
from ullyses_utils.ullyses_config import VERSION
from ullyses_utils.readwrite_yaml import read_config

UTILS_DIR = ullyses_utils.__path__[0]


def read_tss_yaml(targ, yamlfile=None):
    if yamlfile is None:
        yamlfile = os.path.join(UTILS_DIR, f"data/timeseries/{targ}.yaml")
        assert os.path.exists(yamlfile), f"YAML file not found for {targ}, expected {yamlfile}"

    tss_params = read_config(yamlfile)
    return tss_params


def replace_utils_dir(shift_file, utils_dir=UTILS_DIR):
    spl = shift_file.split("/")
    if "$UTILS_DIR" in spl:
        ind = spl.index("$UTILS_DIR")
        spl[ind] = utils_dir
        shift_file = os.path.join(*spl)
    return shift_file


def copy_origdata(datadir, orig_datadir, tss_params):
    """
    Copy the original raw data, to ensure nothing is mistakenly edited.
    Data is copied from orig_datadir to datadir.

    Args:
        datadir (str): Destination directory to copy data into
        orig_datadir (str): Path to the original raw data
    """

    bad_ipppss = tss_params["bad_ipppss"]
    if bad_ipppss is None:
        bad_ipppss = [None]
    
    files = glob.glob(os.path.join(orig_datadir, "l*corrtag*fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*rawtag*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*spt*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*asn*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*x1d.fits"))
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in files:
        ipppss = os.path.basename(item)[:6]
        if ipppss not in bad_ipppss:
            shutil.copy(item, datadir)
    print(f"\nCopied original files to {datadir}")

    return datadir


def calibrate_data(datadir, tss_params, custom_caldir=None): 
    """
    Calibrate data which require special calibration. This will be for
    any exposures which were offset in wavelength

    Args:
        datadir (str): Path to directory with raw timeseries datasets
        custom_caldir (str): Path to output directory

    Returns:
        custom_caldir: Path to output directory with custom data products
    """

    calrequired0 = []
    
    wl_shift_dict = tss_params["wavelength_shift"]
    if wl_shift_dict is not None:
        wl_shift_ipppss = list(wl_shift_dict.keys())
        calrequired0 += wl_shift_ipppss
    else:
        wl_shift_ipppss = [None]
    calrequired = [x[:6] for x in calrequired0]

    bad_ipppss = tss_params["bad_ipppss"]
    if bad_ipppss is None:
        bad_ipppss = [None]

    asns = glob.glob(os.path.join(datadir, "*asn.fits"))
    if custom_caldir is None:
        custom_caldir = os.path.join(datadir, "custom_calibration")
    if not os.path.exists(custom_caldir):
        os.makedirs(custom_caldir)
    for asn in asns:
        ipppss = os.path.basename(asn)[:6]
        if ipppss not in calrequired:
            continue
        if ipppss in bad_ipppss: # This is a bad visit
            continue
        if ipppss in wl_shift_ipppss:
            shift_file0 = wl_shift_dict[ipppss]
            shift_file = replace_utils_dir(shift_file0) 
            calcos.calcos(asn, shift_file=shift_file, outdir=custom_caldir, verbosity=0)
        else:
            calcos.calcos(asn, outdir=custom_caldir, verbosity=0)

    print(f"\nCalibrated data that required a wavelength shift\n\tOutput dir={custom_caldir}")

    return custom_caldir


def copy_caldata(datadir, tss_params, custom_caldir=None):
    """
    Copy the products for each target.

    Args:
        datadir (str): Path to directory with timeseries datasets
        custom_caldir (str): Path to directory with custom calcos outputs

    """
    
    if custom_caldir is None:
        custom_caldir = os.path.join(datadir, "custom_calibration")

    # Copy the custom calibrated output
    corrs = glob.glob(os.path.join(custom_caldir, "*corrtag*"))
    datadir_fuv = os.path.join(datadir, "g160m")
    datadir_nuv = os.path.join(datadir, "g230l")
    datadir_fuvx1d = os.path.join(datadir_fuv, "exp")
    datadir_nuvx1d = os.path.join(datadir_nuv, "exp")
    for d in [datadir_fuv, datadir_nuv, datadir_fuvx1d, datadir_nuvx1d]:
        if not os.path.exists(d):
            os.makedirs(d)
    for item in corrs:
        if fits.getval(item, "opt_elem") == "G160M":
            d = datadir_fuv
        else:
            d = datadir_nuv
#        filename = os.path.basename(item)
#        sptfile = filename.split("_")[0]+"_spt.fits"
#        spt = os.path.join(datadir, sptfile)
    # First copy corrtags, to be turned into splittags
        shutil.copy(item, d)
#        shutil.copy(spt, d)
        
    # Then copy x1ds, which are used as is
        root = os.path.basename(item)[:9]
        x1ds = glob.glob(os.path.join(custom_caldir, root+"*x1d.fits"))
        for x1d in x1ds:
            shutil.copy(x1d, os.path.join(d, "exp"))

    bad_ipppss = tss_params["bad_ipppss"]
    if bad_ipppss is None:
        bad_ipppss = [None]
    
    # Then copy the corrtags that didn't require custom calibration
    orig_corrs = glob.glob(os.path.join(datadir, "*corrtag*"))
    orig_corrfiles = [os.path.basename(x) for x in orig_corrs]
    corrfiles = [os.path.basename(x) for x in corrs]
    for i in range(len(orig_corrs)):
        ipppss = os.path.basename(orig_corrs[i])[:6]
        if orig_corrfiles[i] not in corrfiles and ipppss not in bad_ipppss:
            if fits.getval(orig_corrs[i], "opt_elem") == "G160M":
                d = datadir_fuv
            else:
                d = datadir_nuv
            shutil.copy(orig_corrs[i], d)

    x1ds = glob.glob(os.path.join(custom_caldir, "*x1d.fits"))
    orig_x1ds = glob.glob(os.path.join(datadir, "*x1d.fits"))
    orig_x1dfiles = [os.path.basename(x) for x in orig_x1ds]
    x1dfiles = [os.path.basename(x) for x in x1ds]
    for i in range(len(orig_x1ds)):
        ipppss = os.path.basename(orig_x1ds[i])[:6]
        if orig_x1dfiles[i] not in x1dfiles and ipppss not in bad_ipppss: 
            if fits.getval(orig_x1ds[i], "opt_elem") == "G160M":
                d = datadir_fuv
            else:
                d = datadir_nuv
            shutil.copy(orig_x1ds[i], os.path.join(d, "exp"))

    print(f"\nCopied all corrtags and x1ds to\n\t {datadir}/g160m/ and g230l/\n")
    

def create_splittags(datadir, tss_params):
    """
    Split the input data into discrete time intervals
    by calling splittag_wrapper.

    Args:
        datadir (str): Path to data directory
    """

    bins = tss_params["bins"]
    for grat in ["g160m", "g230l"]:
        for epoch_ippp in bins:
            splitdir = os.path.join(datadir, grat, f"{epoch_ippp}_split")
            time_bin = bins[epoch_ippp][grat]["time"]
            splittag_wrapper.main(indir=splitdir, outdir=splitdir, incr=time_bin, 
                                  numcores=10)
            # The default output directory of the calibrated split corrtags is calcosout
            # We need to move the x1ds from that default directory down a level
            # to the split directory
            calcos_out = glob.glob(os.path.join(splitdir, "calcosout", "*x1d.fits")) 
            for item in calcos_out:
                shutil.move(item, splitdir)


def correct_vignetting(datadir):
    """
    Correct for vignetting on the COS/NUV detector in the G230L/2950 mode.

    Args:
        datadir (str):

    Returns:
        None
    """

    indirs = [os.path.join(datadir, "g230l", "exp"),
              os.path.join(datadir, "g230l", "split")]
    for indir in indirs:
        files = glob.glob(os.path.join(indir, "*x1d.fits"))
        for item in files:
            if fits.getval(item, "cenwave") == 2950:
                root = fits.getval(item, "rootname").lower()
                scale_file = os.path.join(UTILS_DIR, "data/vignette_scaling", f"{root}_scale.txt")
                assert os.path.exists(scale_file), f"No scaling file found for {item}"
                scale = np.loadtxt(scale_file)
                with fits.open(item, mode="update") as hdulist:
                    assert len(scale) == len(hdulist[1].data["flux"][1]),\
                        f"Shape of FITS and scaling factor do not match for {item}"
                    hdulist[1].data["flux"][1] /= scale  # NUVB is 0th index

    print(f'\nApplied scaling factor to G230L/2950 NUVB data in {os.path.join(datadir, "g230l")}') 


def create_timeseries(datadir, tss_outdir, targ, tss_params):
    """
    Creates the timeseries high level science products for ULLYSES monitoring targets
    from the custom-calibrated and split data products.

    Args:
        datadir (str): Path to input data
        tss_outdir (str): Path to time series data product output directory
        targ (str): Target name
    """

    bins = tss_params["bins"]
    for grat in ["g160m", "g230l"]:
        epoch_ippp = list(bins.keys())[0]
        wl_bin = bins[epoch_ippp][grat]["wave"]
        min_exptime = bins[epoch_ippp][grat]["min_exptime"]
        if not os.path.exists(tss_outdir):
            os.makedirs(tss_outdir)
        # First create the exposure level time-series spectra
        indir = os.path.join(datadir, grat, "exp")
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{VERSION.lower()}_tss.fits")
        timeseries.process_files(grat.upper(), outfile, indir, overwrite=True) 
        
        indir = os.path.join(datadir, grat, "split")
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{VERSION.lower()}_split-tss.fits")
        timeseries.process_files(grating=grat.upper(), outfile=outfile, indir=indir, 
                                 wavelength_binning=wl_bin, min_exptime=min_exptime, 
                                 overwrite=True)

    return tss_outdir


def move_input_epoch_data(datadir, tss_params):
    """
    Each monitoring star includes two epochs of data. The corrtags for each
    epoch need to be kept in separate directories because each epoch may have
    different time sampling. The code to create splittags takes as input a
    single directory, so to create splittags of different time intervals as a
    function of epoch, each epoch's corrtags must be in separate directories.
    """

    bins = tss_params["bins"]
    for grat in ["g160m", "g230l"]:
        singledir = os.path.join(datadir, grat)
        for epoch_ippp in bins:
            splitdir = os.path.join(singledir, f"{epoch_ippp}_split")
            if not os.path.exists(splitdir):
                os.makedirs(splitdir)
            # First move corrtags
            corrtags = glob.glob(os.path.join(singledir, f"{epoch_ippp}*corrtag*fits"))
            for item in corrtags:
                shutil.move(item, splitdir)
            # Now copy SPT files
            for corr in corrtags:
                root = os.path.basename(corr)[:9]
                spt = root+"_spt.fits"
                copied = os.path.join(splitdir, spt)
                if not os.path.exists(copied):
                    sptfile = os.path.join(datadir, spt)
                    shutil.copy(sptfile, splitdir)
    print(f"\nMoved corrtags to be split")


def move_output_epoch_data(datadir, tss_params):
    """
    Once splittags have been made and calibrated for each epoch (where each 
    epoch may have different split time intervals), the resulting x1ds
    must be copied to a single directory to be made into TSS products. 
    The timeseries creation code requires all input files to be stored in a 
    single directory.
    """
    
    bins = tss_params["bins"]
    for grat in ["g160m", "g230l"]:
        single_splitdir = os.path.join(datadir, grat, "split")
        if not os.path.exists(single_splitdir):
            os.makedirs(single_splitdir)
        for epoch_ippp in bins:
            splitdir = os.path.join(datadir, grat, f"{epoch_ippp}_split")
            x1ds = glob.glob(os.path.join(splitdir, "*x1d.fits"))
            for item in x1ds:
                shutil.move(item, single_splitdir)
    print(f"\nMoved x1ds to a single directory to be made into TSS")


def main(datadir, orig_datadir, tss_outdir, targ, yamlfile=None, custom_caldir=None):
    """
    Performs the timeseries data product calibration and creation.

    Args:
        datadir (str):
        orig_datadir (str):
        tss_outdir (str):
        targ (str):
        custom_caldir (str):
    """

    tss_params = read_tss_yaml(targ, yamlfile)
    copy_origdata(datadir, orig_datadir, tss_params)
    calibrate_data(datadir, tss_params, custom_caldir)
    copy_caldata(datadir, tss_params, custom_caldir)
    move_input_epoch_data(datadir, tss_params)
    create_splittags(datadir, tss_params)
    move_output_epoch_data(datadir, tss_params)
    correct_vignetting(datadir)
    create_timeseries(datadir, tss_outdir, targ, tss_params)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--orig",
                        help="Path to all data to be made into a timeseries spectrum")
    parser.add_argument("--copydir",
                        help="Path to copy intermediate products to")
    parser.add_argument("--hlspdir", dest="tss_outdir",
                        help="Directory to write final timeseries HLSPs")
    parser.add_argument("-t", "--targ",
                        help="Name of target")
    parser.add_argument("-y", "--yaml",
                        help="Name of YAML configuration file")
    parser.add_argument("--custom_caldir",
                        help="Name of directory to write splittag products, by default uses current runtime")
    args = parser.parse_args()

    main(args.copydir, args.orig, args.tss_outdir, args.targ, args.yaml, args.custom_caldir)
