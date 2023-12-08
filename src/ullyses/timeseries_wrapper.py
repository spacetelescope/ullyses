"""
Wrapper to make timeseries products, both subexposure and exposure level.
If exposure level -> the input should be the 1D spectra
If subexposure level -> the input should be raw + spt + asn files, so that the files can be split
and recalibrated to make subexposure x1ds.
"""

import traceback
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
from ullyses_utils.readwrite_yaml import read_config
from . import __release__

UTILS_DIR = ullyses_utils.__path__[0]
RED = "\033[1;31m"
RESET = "\033[0;0m"


def remove_identifier(all_ids, ipppss, ipppssoot):
    for identifier in [ipppss, ipppssoot]:
        try:
            all_ids.remove(identifier)
        except:
            pass



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
    assert "good_files" in tss_params, "Good files must be specified for all config files!"
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


def read_tss_yaml(targ, yamlfile=None, instrument="cos"):
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
        yamlfile = os.path.join(UTILS_DIR, f"data/timeseries/{targ}_{instrument.lower()}.yaml")
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


def copy_exp_origdata(datadir, orig_datadir, tss_params):
    """
    For exposure-level data.
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

    should_be_copied = tss_params["good_ipppssoot"] + tss_params["good_ipppss"]
    # Covers both STIS & COS 
    allfiles = glob.glob(os.path.join(orig_datadir, "*x1d.fits"))
    allfiles += glob.glob(os.path.join(orig_datadir, "*sx1.fits"))
    allfiles += glob.glob(os.path.join(orig_datadir, "*x1f.fits"))
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    for item in allfiles:
        basen = os.path.basename(item)
        ipppssoot = basen[:9]
        ipppss = basen[:6]
        if ipppss not in tss_params["good_ipppss"]:
            if ipppssoot not in tss_params["good_ipppssoot"]:
                continue
        if ipppss in tss_params["bad_ipppss"] or ipppssoot in tss_params["bad_ipppssoot"]:
            continue
        remove_identifier(should_be_copied, ipppss, ipppssoot)
        shutil.copy(item, datadir)
    print(f"\nCopied original files to {datadir}")
    if len(should_be_copied) != 0:
        print(f"{RED}WARNING: Not all good files were copied! Missing files: {should_be_copied}{RESET}")

    return datadir


def copy_subexp_origdata(datadir, orig_datadir, tss_params):
    """
    For stars with sub-exposure data.
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

    do_copy = []
    do_copy_ipppss = []
    should_be_copied = tss_params["good_ipppssoot"] + tss_params["good_ipppss"]

    bins_ippp = tss_params["bins"]
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in raws:
        basen = os.path.basename(item)
        ipppssoot = basen[:9]
        ipppss = basen[:6]
        ippp = basen[:4]
        if ipppss not in tss_params["good_ipppss"]:
            if ipppssoot not in tss_params["good_ipppssoot"]:
                continue
        if ipppss in tss_params["bad_ipppss"] or ipppssoot in tss_params["bad_ipppssoot"]:
            continue
        ins = fits.getval(item, "instrume")
        if ins.lower() != tss_params["instrument"]:
            remove_identifier(should_be_copied, ipppss, ipppssoot)
            continue
        grating = fits.getval(item, "opt_elem")
        if grating.lower() not in tss_params["gratings"]:
            remove_identifier(should_be_copied, ipppss, ipppssoot)
            continue
        if ippp not in bins_ippp:
            remove_identifier(should_be_copied, ipppss, ipppssoot)
            continue
        do_copy_ipppss.append(ipppss)
        do_copy += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*rawtag*.fits"))
        do_copy += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*corrtag*.fits"))
        do_copy += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*spt*.fits"))
        do_copy += glob.glob(os.path.join(orig_datadir, f"{ipppssoot}*x1d*.fits"))
        remove_identifier(should_be_copied, ipppss, ipppssoot)

    do_copy = list(set(do_copy))
    do_copy_ipppss = list(set(do_copy_ipppss))

    for item in do_copy:
        shutil.copy(item, datadir)
    asns = glob.glob(os.path.join(orig_datadir, "l*asn*.fits"))
    for item in asns:
        basen = os.path.basename(item)
        ipppss = basen[:6]
        if ipppss in do_copy_ipppss:
            shutil.copy(item, datadir)

    print(f"\nCopied original files to {datadir}")
    if len(should_be_copied) != 0:
        print(f"{RED}WARNING: Not all good files were copied! Missing files: {should_be_copied}{RESET}")
    
    return datadir


def calibrate_cos_data(datadir, tss_params, custom_caldir=None, overwrite=True): 
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

    asns = glob.glob(os.path.join(datadir, "l*asn.fits"))
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

        # Delete existing data if overwrite is True
        if overwrite is True:
            asndata = fits.getdata(asn)
            members = [x.lower() for x in asndata["memname"]]
            for member in members:
                outfiles = glob.glob(os.path.join(custom_caldir, member+"*"))
                if len(outfiles) > 0:
                    print(f"Overwrite is True, removing existing products for {member}...")
                    for outfile in outfiles:
                        os.remove(outfile)

        if ipppss in wl_shift_ipppss:
            shift_file0 = wl_shift_dict[ipppss]
            shift_file = replace_utils_dir(shift_file0) 
            calcos.calcos(asn, shift_file=shift_file, outdir=custom_caldir, verbosity=0)
        else:
            calcos.calcos(asn, outdir=custom_caldir, verbosity=0)

    print(f"\nCalibrated data that required a wavelength shift\n\tOutput dir={custom_caldir}")

    return custom_caldir


def copy_subexp_caldata(datadir, tss_params, custom_caldir=None):
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
        elif fits.getval(item, "opt_elem") == "G230L" and fits.getval(item, "instrume") == "COS":
            d = datadir_g230l
        else:
            continue
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
            elif fits.getval(orig_corrs[i], "opt_elem") == "G230L" and fits.getval(orig_corrs[i], "instrume") == "COS":
                d = datadir_g230l
            else:
                continue
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
            elif fits.getval(orig_x1ds[i], "opt_elem") == "G230L" and fits.getval(orig_x1ds[i], "instrume") == "COS":
                d = datadir_g230l
            else:
                continue
            shutil.copy(orig_x1ds[i], os.path.join(d, "exp"))

    print(f"\nCopied all corrtags and x1ds to\n\t {datadir}/g160m/ and g230l/\n")
    

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
            if grat not in tss_params["bins"][epoch_ippp]:
                continue
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
                try:
                    scl_done = fits.getval(item, "scl_done")
                except KeyError:
                    scl_done = "False"
                if scl_done == "True":
                    continue
                root = fits.getval(item, "rootname").lower()
                scale_file = os.path.join(UTILS_DIR, "data/vignette_scaling", f"{root}_scale.txt")
                assert os.path.exists(scale_file), f"No scaling file found for {item}"
                scale = np.loadtxt(scale_file)
                with fits.open(item, mode="update") as hdulist:
                    assert len(scale) == len(hdulist[1].data["flux"][1]),\
                        f"Shape of FITS and scaling factor do not match for {item}"
                    hdulist[1].data["flux"][1] /= scale  # NUVB is 1st index
                    hdulist[0].header["SCL_DONE"] = "True" 

    print(f'\nApplied scaling factor to G230L/2950 NUVB data in directories: \n{indirs}\n') 


def create_exp_timeseries(datadir, tss_outdir, targ, tss_params, min_exptime=0.1):
    """
    Creates the timeseries high level science products for targets with exposure-level TSS
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
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_{ins}_{targ}_{grat}_{__release__.lower()}_tss.fits")
        try:
            timeseries.process_files(grat.upper(), outfile, datadir, overwrite=True, 
                                 ins=ins.upper(), min_exptime=min_exptime)
        except Exception:
            print(f"{RED}WARNING: EXCEPTION={RESET}")
            without = glob.glob(os.path.join(datadir, "*without*fits"))
            for newname in without:
                origname = newname.replace("_without.fits", "_x1d.fits")
                os.rename(newname, origname)
            print(traceback.format_exc())

    return tss_outdir


def create_subexp_timeseries(datadir, tss_outdir, targ, tss_params):
    """
    Creates the timeseries high level science products for with sub-exposure level TSS
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
        if grat not in ["g230l", "g160m"]:
            continue
        if not os.path.exists(tss_outdir):
            os.makedirs(tss_outdir)
        # First create the exposure level time-series spectra
        indir = os.path.join(datadir, grat, "exp")
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{__release__.lower()}_tss.fits")
        try:
            timeseries.process_files(grat.upper(), outfile, indir, overwrite=True, min_exptime=1) 
        except Exception:
            print(f"{RED}WARNING: EXCEPTION={RESET}")
            without = glob.glob(os.path.join(datadir, "*without*fits"))
            for newname in without:
                origname = newname.replace("_without.fits", "_x1d.fits")
                os.rename(newname, origname)
            print(traceback.format_exc())
       
        # Now make sub-exposure level time series
        epoch_ippp = list(bins.keys())[0]
        if grat not in tss_params["bins"][epoch_ippp]:
            continue
        wl_bin = bins[epoch_ippp][grat]["wave"]
        min_exptime = bins[epoch_ippp][grat]["min_exptime"]
        indir = os.path.join(datadir, grat, "split")
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{__release__.lower()}_split-tss.fits")
        try:
            timeseries.process_files(grating=grat.upper(), outfile=outfile, indir=indir, 
                                 wavelength_binning=wl_bin, min_exptime=min_exptime, 
                                 overwrite=True)
        except Exception:
            print(f"{RED}WARNING: EXCEPTION={RESET}")
            without = glob.glob(os.path.join(datadir, "*without*fits"))
            for newname in without:
                origname = newname.replace("_without.fits", "_x1d.fits")
                os.rename(newname, origname)
            print(traceback.format_exc())

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
        if grat not in ["g160m", "g230l"]:
            continue
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
        print(f"\nMoved corrtags to be split: {splitdir}")


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
        if grat not in ["g160m", "g230l"]:
            continue
        single_splitdir = os.path.join(datadir, grat, "split")
        if not os.path.exists(single_splitdir):
            os.makedirs(single_splitdir)
        for epoch_ippp in bins:
            splitdir = os.path.join(datadir, grat, f"{epoch_ippp}_split")
            x1ds = glob.glob(os.path.join(splitdir, "*x1d.fits"))
            for item in x1ds:
                shutil.move(item, single_splitdir)
        print(f"\nMoved x1ds to a single directory to be made into TSS: {single_splitdir}")


def exp_star(datadir, orig_datadir, tss_outdir, targ, yamlfile=None, custom_caldir=None,
                       min_exptime=0.1, instrument="cos"):
    """
    Args:
        datadir (str): The path the original data will be copied to.
        orig_datadir (str): Path to the original raw data
        tss_outdir (str): Path to time series data product output directory
        targ (str): Name of target. Must be the ULLYSES DP target name!!!
        custom_caldir (str): Path to directory with custom calcos outputs
    """
    if not os.path.exists(tss_outdir):
        os.makedirs(tss_outdir)
    tss_params = read_tss_yaml(targ, yamlfile, instrument)
    tss_params = get_goodbad_exposures(tss_params)
    copy_exp_origdata(datadir, orig_datadir, tss_params)
    create_exp_timeseries(datadir, tss_outdir, targ, tss_params, min_exptime)


def subexp_star(datadir, orig_datadir, tss_outdir, targ, yamlfile=None, custom_caldir=None, 
                    instrument="cos"):
    """
    Main wrapper function to create TSS for a star with sub-exposure TSS.
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
    tss_params = read_tss_yaml(targ, yamlfile, instrument)
    tss_params = get_goodbad_exposures(tss_params)
    copy_subexp_origdata(datadir, orig_datadir, tss_params)
    calibrate_cos_data(datadir, tss_params, custom_caldir)
    copy_subexp_caldata(datadir, tss_params, custom_caldir)
    move_input_epoch_data(datadir, tss_params)
    create_splittags(datadir, tss_params)
    move_output_epoch_data(datadir, tss_params)
    if "g230l" in tss_params["gratings"]:
        correct_vignetting(datadir)
    create_subexp_timeseries(datadir, tss_outdir, targ, tss_params)


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
