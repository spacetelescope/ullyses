"""
Wrapper to make timeseries products, both subexposure and exposure level.
If exposure level -> the input should be the 1D spectra
If subexposure level -> the input should be raw + spt + asn files, so that the files can be split
and recalibrated to make subexposure x1ds.
"""

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


def get_goodbad_exposures(tss_params):
    """Determine if there are any bad exposures to avoid.

    Args:
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.

    Returns:
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use. This has been
            updated to explicitly list bad IPPPSS and IPPPSSOOT separately.
    """
    for filetype in ["bad", "good"]:
        filetype_files = tss_params[f"{filetype}_files"]
        filetype_ipppss = []
        filetype_ipppssoot = []
        if filetype_files is None:
            filetype_ipppss.append(None)
            filetype_ipppssoot.append(None)
        else:
            for item in filetype_files:
                if len(item) == 9:
                    filetype_ipppssoot.append(item)
                elif len(item) == 7 and "*" in item:
                    ipppss = item.strip("*")
                    filetype_ipppss.append(ipppss)
                elif len(item) == 6:
                    filetype_ipppss.append(item)
        tss_params[f"{filetype}_ipppss"] = filetype_ipppss
        tss_params[f"{filetype}_ipppssoot"] = filetype_ipppssoot
    
    return tss_params


def read_tss_yaml(targ, yamlfile=None):
    """ Open the YAML file for the target to read TSS parameters.

    Args:
        targ (str): Name of target. Must be the ULLYSES DP target name!!!
        yamlfile (str or None): (Optional) The name of the yamlfile to use for
            the supplied target. If none is specified, it will be looked up
            on the fly from the ullyses_utils repository.

    Returns:
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    """

    if yamlfile is None:
        yamlfile = os.path.join(UTILS_DIR, f"data/timeseries/{targ}.yaml")
        assert os.path.exists(yamlfile), f"YAML file not found for {targ}, expected {yamlfile}"

    tss_params = read_config(yamlfile)
    return tss_params


def replace_utils_dir(shift_file, utils_dir=UTILS_DIR):
    """Replace the $UTILS_DIR variable with the absolute path.

    Args:
        shift_file (str): The path to the wavelength shift file which includes
            a $UTILS_DIR reference.
        utils_dir (str): The absolute path to the ullyses_utils repository.
    Returns:
        shift_file (str): Updated shift_file variable with the absolute utils
            path.
    """

    spl = shift_file.split("/")
    if "$UTILS_DIR" in spl:
        ind = spl.index("$UTILS_DIR")
        spl[ind] = utils_dir
        shift_file = os.path.join(*spl)
    return shift_file


def copy_serendipitous_origdata(datadir, orig_datadir, tss_params):
    """
    NOT CURRENTLY IN USE.
    For serendipitous stars.
    Copy the original raw data, to ensure nothing is mistakenly edited.
    Data is copied from orig_datadir to datadir.

    Args:
        datadir (str): Destination directory to copy data into
        orig_datadir (str): Path to the original raw data
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    Returns:
        datadir (str): The path where data were copied to. Sometimes this
            directory is created on the fly. 
    """

    # Covers both STIS & 
    files = glob.glob(os.path.join(orig_datadir, "*corrtag*fits"))
    files += glob.glob(os.path.join(orig_datadir, "*rawtag*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "*spt*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "*asn*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "*x1d.fits"))
    files += glob.glob(os.path.join(orig_datadir, "*raw.fits"))
    files += glob.glob(os.path.join(orig_datadir, "*sx1.fits"))
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in files:
        basen = os.path.basename(item)
        ipppssoot = basen[:9]
        ipppss = basen[:6]
        if ipppss in tss_params["bad_ipppss"] or ipppssoot in tss_params["bad_ipppssoot"]:
            continue
        shutil.copy(item, datadir)
    print(f"\nCopied original files to {datadir}")

    return datadir


def copy_monitoring_origdata(datadir, orig_datadir, tss_params):
    """
    For monitoring stars.
    Copy the original raw data, to ensure nothing is mistakenly edited.
    Data is copied from orig_datadir to datadir.

    Args:
        datadir (str): Destination directory to copy data into
        orig_datadir (str): Path to the original raw data
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    Returns:
        datadir (str): The path where data were copied to. Sometimes this
            directory is created on the fly. 
    """

    raws = glob.glob(os.path.join(orig_datadir, "l*rawtag*.fits"))

    good_ipppss = []
    good_ipppssoot = []

    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in raws:
        basen = os.path.basename(item)
        ipppssoot = basen[:9]
        ipppss = basen[:6]
        if ipppss in tss_params["bad_ipppss"] or ipppssoot in tss_params["bad_ipppssoot"]:
            continue
        ins = fits.getval(item, "instrume")
        if ins.lower() != tss_params["instrument"]:
            continue
        grating = fits.getval(item, "opt_elem")
        if grating.lower() not in tss_params["gratings"]:
            continue
        good_ipppss.append(ipppss)
        good_ipppssoot += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*rawtag*.fits"))
        good_ipppssoot += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*corrtag*.fits"))
        good_ipppssoot += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*spt*.fits"))
        good_ipppssoot += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*x1d*.fits"))

    good_ipppssoot = list(set(good_ipppssoot))
    good_ipppss = list(set(good_ipppss))

    for item in good_ipppssoot:
        shutil.copy(item, datadir)
    asns = glob.glob(os.path.join(orig_datadir, "l*asn*.fits"))
    for item in asns:
        basen = os.path.basename(item)
        ipppss = basen[:6]
        if ipppss in good_ipppss:
            shutil.copy(item, datadir)

    print(f"\nCopied original files to {datadir}")

    return datadir


def calibrate_cos_data(datadir, tss_params, custom_caldir=None): 
    """
    Calibrate COS data which require special calibration. Typically for
    any exposures which were offset in wavelength.

    Args:
        datadir (str): Path to directory with raw timeseries datasets
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
        custom_caldir (str): Path to output directory

    Returns:
        custom_caldir: Path to output directory with custom data products.
            Sometimes this directory is created on the fly.
    """

    calrequired0 = []
    
    wl_shift_dict = tss_params["wavelength_shift"]
    if wl_shift_dict is not None:
        wl_shift_ipppss = list(wl_shift_dict.keys())
        calrequired0 += wl_shift_ipppss
    else:
        wl_shift_ipppss = [None]
    calrequired = [x[:6] for x in calrequired0]

    asns = glob.glob(os.path.join(datadir, "*asn.fits"))
    if custom_caldir is None:
        custom_caldir = os.path.join(datadir, "custom_calibration")
    if not os.path.exists(custom_caldir):
        os.makedirs(custom_caldir)
    for asn in asns:
        ipppss = os.path.basename(asn)[:6]
        if ipppss not in calrequired:
            continue
        if ipppss in tss_params["bad_ipppss"]: # This is a bad visit
            continue
        if ipppss in wl_shift_ipppss:
            shift_file0 = wl_shift_dict[ipppss]
            shift_file = replace_utils_dir(shift_file0) 
            calcos.calcos(asn, shift_file=shift_file, outdir=custom_caldir, verbosity=0)
        else:
            calcos.calcos(asn, outdir=custom_caldir, verbosity=0)

    print(f"\nCalibrated data that required a wavelength shift\n\tOutput dir={custom_caldir}")

    return custom_caldir


def copy_monitoring_caldata(datadir, tss_params, custom_caldir=None):
    """
    Copy the calibrated products for each target (corrtags and x1ds).

    Args:
        datadir (str): Path to directory with timeseries datasets
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
        custom_caldir (str): Path to directory with custom calcos outputs

    """
    
    if custom_caldir is None:
        custom_caldir = os.path.join(datadir, "custom_calibration")

    # Copy any custom calibrated output
    corrs = glob.glob(os.path.join(custom_caldir, "*corrtag*"))
    datadir_g160m = os.path.join(datadir, "g160m")
    datadir_g230l = os.path.join(datadir, "g230l")
    datadir_g160mx1d = os.path.join(datadir_g160m, "exp")
    datadir_g230lx1d = os.path.join(datadir_g230l, "exp")
    for d in [datadir_g160m, datadir_g230l, datadir_g160mx1d, datadir_g230lx1d]:
        if not os.path.exists(d):
            os.makedirs(d)
    for item in corrs:
        if fits.getval(item, "opt_elem") == "G160M":
            d = datadir_g160m
        else:
            d = datadir_g230l
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

    # Then copy the corrtags that didn't require custom calibration
    orig_corrs = glob.glob(os.path.join(datadir, "*corrtag*"))
    orig_corrfiles = [os.path.basename(x) for x in orig_corrs]
    corrfiles = [os.path.basename(x) for x in corrs]
    for i in range(len(orig_corrs)):
        ipppss = os.path.basename(orig_corrs[i])[:6]
        ipppssoot = os.path.basename(orig_corrs[i])[:9]
        if orig_corrfiles[i] not in corrfiles and ipppss not in tss_params["bad_ipppss"] and ipppssoot not in tss_params["bad_ipppssoot"]:
            if fits.getval(orig_corrs[i], "opt_elem") == "G160M":
                d = datadir_g160m
            else:
                d = datadir_g230l
            shutil.copy(orig_corrs[i], d)

    x1ds = glob.glob(os.path.join(custom_caldir, "*x1d.fits"))
    orig_x1ds = glob.glob(os.path.join(datadir, "*x1d.fits"))
    orig_x1dfiles = [os.path.basename(x) for x in orig_x1ds]
    x1dfiles = [os.path.basename(x) for x in x1ds]
    for i in range(len(orig_x1ds)):
        ipppss = os.path.basename(orig_x1ds[i])[:6]
        ipppssoot = os.path.basename(orig_x1ds[i])[:9]
        if orig_x1dfiles[i] not in x1dfiles and ipppss not in tss_params["bad_ipppss"] and ipppssoot not in tss_params["bad_ipppssoot"]:
            if fits.getval(orig_x1ds[i], "opt_elem") == "G160M":
                d = datadir_g160m
            else:
                d = datadir_g230l
            shutil.copy(orig_x1ds[i], os.path.join(d, "exp"))

    print(f"\nCopied all corrtags and x1ds to\n\t {datadir}g160m/ and g230l/\n")
    

def create_splittags(datadir, tss_params):
    """
    Split the input data into discrete time intervals
    by calling splittag_wrapper.

    Args:
        datadir (str): Path to data directory
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    """

    bins = tss_params["bins"]
    for grat in tss_params["gratings"]:
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
        datadir (str): Path to data directory.

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
                    hdulist[1].data["flux"][1] /= scale  # NUVB is 1st index

    print(f'\nApplied scaling factor to G230L/2950 NUVB data in {os.path.join(datadir, "g230l")}\n') 


def create_serendipitous_timeseries(datadir, tss_outdir, targ, tss_params, min_exptime=1):
    """
    Creates the timeseries high level science products for ULLYSES monitoring targets
    from the custom-calibrated and split data products.

    Args:
        datadir (str): Path to input data
        tss_outdir (str): Path to time series data product output directory
        targ (str): Target name
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    Returns:
        tss_outdir (str): Path to input data, sometimes this is created on the fly.
    """

    ins = tss_params["instrument"]
    for grat in tss_params["gratings"]: 
        # Create the exposure level time-series spectra
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_{ins}_{targ}_{grat}_{VERSION.lower()}_tss.fits")
        timeseries.process_files(grat.upper(), outfile, datadir, overwrite=True, 
                                 ins=ins.upper(), min_exptime=min_exptime) 
        
    return tss_outdir



def create_monitoring_timeseries(datadir, tss_outdir, targ, tss_params):
    """
    Creates the timeseries high level science products for ULLYSES monitoring targets
    from the custom-calibrated and split data products.

    Args:
        datadir (str): Path to input data
        tss_outdir (str): Path to time series data product output directory
        targ (str): Target name
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    Returns:
        tss_outdir (str): Path to input data, sometimes this is created on the fly.
    """

    bins = tss_params["bins"]
    for grat in tss_params["gratings"]:
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
    Move data from different epochs into different directories.
    Each monitoring star includes two epochs of data. The corrtags for each
    epoch need to be kept in separate directories because each epoch may have
    different time sampling. The code to create splittags takes as input a
    single directory, so to create splittags of different time intervals as a
    function of epoch, each epoch's corrtags must be in separate directories.

    Args:
        datadir (str): Path to input data
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    """

    bins = tss_params["bins"]
    for grat in tss_params["gratings"]:
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

    Args:
        datadir (str): Path to input data
        tss_params (dict): A dictionary where keys are the TSS parameters and
            values are the parameter values, e.g. files to use.
    """
    
    bins = tss_params["bins"]
    for grat in tss_params["gratings"]:
        single_splitdir = os.path.join(datadir, grat, "split")
        if not os.path.exists(single_splitdir):
            os.makedirs(single_splitdir)
        for epoch_ippp in bins:
            splitdir = os.path.join(datadir, grat, f"{epoch_ippp}_split")
            x1ds = glob.glob(os.path.join(splitdir, "*x1d.fits"))
            for item in x1ds:
                shutil.move(item, single_splitdir)
    print(f"\nMoved x1ds to a single directory to be made into TSS")


def serendipitous_star(datadir, tss_outdir, targ, yamlfile=None, custom_caldir=None,
                       min_exptime=1):
    if not os.path.exists(tss_outdir):
        os.makedirs(tss_outdir)
    tss_params = read_tss_yaml(targ, yamlfile)
    tss_params = get_goodbad_exposures(tss_params)
    # Not yet implemented because it's not needed so far.
    #copy_serendipitous_origdata(datadir, orig_datadir, tss_params)
    #calibrate_cos_data(datadir, tss_params, custom_caldir)
    #copy_monitoring_caldata(datadir, tss_params, custom_caldir)
    create_serendipitous_timeseries(datadir, tss_outdir, targ, tss_params, min_exptime)


def monitoring_star(datadir, orig_datadir, tss_outdir, targ, yamlfile=None, custom_caldir=None):
    """
    Main wrapper function to create TSS for a monitoring star.
    Performs the timeseries data product calibration and creation.

    Args:
        datadir (str): The path the original data will be copied to.
        orig_datadir (str): Path to the original raw data
        tss_outdir (str): Path to time series data product output directory
        targ (str): Name of target. Must be the ULLYSES DP target name!!!
        custom_caldir (str): Path to directory with custom calcos outputs
    """

    if not os.path.exists(tss_outdir):
        os.makedirs(tss_outdir)
    tss_params = read_tss_yaml(targ, yamlfile)
    tss_params = get_goodbad_exposures(tss_params)
    copy_monitoring_origdata(datadir, orig_datadir, tss_params)
    calibrate_cos_data(datadir, tss_params, custom_caldir)
    copy_monitoring_caldata(datadir, tss_params, custom_caldir)
    move_input_epoch_data(datadir, tss_params)
    create_splittags(datadir, tss_params)
    move_output_epoch_data(datadir, tss_params)
    if "g230l" in tss_params["gratings"]:
        correct_vignetting(datadir)
    create_monitoring_timeseries(datadir, tss_outdir, targ, tss_params)


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
