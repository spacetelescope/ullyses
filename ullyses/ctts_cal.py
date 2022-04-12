import argparse
import numpy as np
from astropy.io import fits
import glob
import os
import shutil
import calcos

import splittag_wrapper
import timeseries
import ullyses_utils
from ullyses_utils.ullyses_config import VERSION

utils_dir = ullyses_utils.__path__[0]

BINS = {"v-tw-hya": {"g160m": {"time": 30, "wave": 3, "min_exptime": 20}, # exptime = 300
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}, # exptime = 30
        "v-ru-lup": {"g160m": {"time": 30, "wave": 6, "min_exptime": 20}, # exptime = 220
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}, # exptime = 30
        "v-bp-tau": {"g160m": {"time": 30, "wave": 6, "min_exptime": 20}, # exptime = 128
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}, # exptime = 196
        "v-gm-aur": {"g160m": {"time": 90, "wave": 6, "min_exptime": 50}, # exptime = 186
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}} # exptime = 184

BAD_IPPPSS = ["le9d1k",  # TW Hydra
              "leit1d", "leitad", "leit1l",  # RU Lup
              "lek71f"]  # GM Aur

G230L_DISPTAB = os.path.join(utils_dir, 'data/ref_files/ullyses_cos_nuv_disp.fits')

WL_SHIFT = {'le9d1c': os.path.join(utils_dir, "data/cos_shifts/twhya_shifts.txt"), 
            'le9d1g': os.path.join(utils_dir, "data/cos_shifts/twhya_shifts.txt")} 
# TW Hydra exposures
#['le9d1cdeq',
# 'le9d1cdeq',
# 'le9d1cdgq',
# 'le9d1cdgq',
# 'le9d1cdiq',
# 'le9d1cdiq',
# 'le9d1cdkq',
# 'le9d1cdkq',
# 'le9d1cdmq',
# 'le9d1cdmq',
# 'le9d1cdoq',
# 'le9d1cdoq',
# 'le9d1gw7q']


def copy_origdata(datadir, orig_datadir):
    """
    Copy the original raw data, to ensure nothing is mistakenly edited.
    Data is copied from orig_datadir to datadir.

    Args:
        datadir (str): Destination directory to copy data into
        orig_datadir (str): Path to the original raw data
    """

    files = glob.glob(os.path.join(orig_datadir, "l*corrtag*fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*rawtag*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*spt*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*asn*.fits"))
    files += glob.glob(os.path.join(orig_datadir, "l*x1d.fits"))
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in files:
        shutil.copy(item, datadir)
    print(f"Copied original files to {datadir}")

    return datadir


def calibrate_data(datadir, custom_caldir=None, g230l_disptab=G230L_DISPTAB):
    """
    Calibrate data which require special calibration. This will be for
    1) all G230L data which need to be calibrated with a custom COS/NUV DISPTAB
    2) some TW Hydra exposures which were offset in wavelength

    Args:
        datadir (str): Path to directory with raw timeseries datasets
        custom_caldir (str): Path to output directory
        g230l_disptab (str): Path to custom COS NUV DISPTAB

    Returns:
        custom_caldir: Path to output directory with custom data products
    """

    calrequired0 = []
    raws = glob.glob(os.path.join(datadir, "l*rawtag*fits"))
    for raw in raws:
        if fits.getval(raw, "opt_elem") == "G230L":
            with fits.open(raw, mode="update") as hdulist:
                hdr0 = hdulist[0].header
                hdr0.set("DISPTAB", g230l_disptab)
            root = os.path.basename(raw)[:9]
            calrequired0.append(root)
    wl_shift_files = list(WL_SHIFT.keys())
    calrequired0 += wl_shift_files
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
        if ipppss in BAD_IPPPSS: # This is a bad visit
            continue
        if ipppss in WL_SHIFT: 
            calcos.calcos(asn, shift_file=WL_SHIFT[ipppss], outdir=custom_caldir, verbosity=0)
        else:
            calcos.calcos(asn, outdir=custom_caldir, verbosity=0)

    print(f"\n Calibrated data that required a shift or custom DISPTAB\n, output dir={custom_caldir}")

    return custom_caldir


def copy_caldata(datadir, custom_caldir=None):
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

    # Then copy the corrtags that didn't require custom calibration
    orig_corrs = glob.glob(os.path.join(datadir, "*corrtag*"))
    orig_corrfiles = [os.path.basename(x) for x in orig_corrs]
    corrfiles = [os.path.basename(x) for x in corrs]
    for i in range(len(orig_corrs)):
        ipppss = os.path.basename(orig_corrs[i])[:6]
        if orig_corrfiles[i] not in corrfiles and ipppss not in BAD_IPPPSS:
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
        if orig_x1dfiles[i] not in x1dfiles and ipppss not in BAD_IPPPSS: 
            if fits.getval(orig_corrs[i], "opt_elem") == "G160M":
                d = datadir_fuv
            else:
                d = datadir_nuv
            shutil.copy(orig_corrs[i], os.path.join(d, "exp"))

    print(f"\nCopied all corrtags and x1ds to {datadir}/g160m/ and g230l/\n")
    

def create_splittags(datadir, targ):
    """
    Split the input data into several different time and wavlength bins
    by calling splittag_wrapper.

    Args:
        datadir (str): Path to data directory
        targ (str): Target name
    """

    for grat in ["g160m", "g230l"]:
        indir = os.path.join(datadir, grat)
        splitdir = os.path.join(indir, "split")
        if not os.path.exists(splitdir):
            os.makedirs(splitdir)
        # First copy SPT files
        infiles = glob.glob(os.path.join(indir, "*corrtag*fits"))
        for infile in infiles:
            root = os.path.basename(infile)[:9]
            spt = root+"_spt.fits"
            copied = os.path.join(splitdir, spt)
            if not os.path.exists(copied):
                sptfile = os.path.join(datadir, spt)
                shutil.copy(sptfile, splitdir)

        t = BINS[targ][grat]["time"]
        splittag_wrapper.main(datadir, splitdir, incr=t, numcores=10)
        # The default output directory of the calibrated split corrtags is calcosout
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
                scale_file = os.path.join(utils_dir, "data/vignette_scaling", f"{root}_scale.txt")
                assert os.path.exists(scale_file), f"No scaling file found for {item}"
                scale = np.loadtxt(scale_file)
                with fits.open(item, mode="update") as hdulist:
                    assert len(scale) == len(hdulist[1].data["flux"][1]),\
                        f"Shape of FITS and scaling factor do not match for {item}"
                    hdulist[1].data["flux"][1] /= scale  # NUVB is 0th index

    print(f'Applied scaling factor to G230L/2950 NUVB data in {os.path.join(datadir, "g230l")}') 


def create_timeseries(datadir, tss_outdir, targ):
    """
    Creates the timeseries high level science products for ULLYSES monitoring targets
    from the custom-calibrated and split data products.

    Args:
        datadir (str): Path to input data
        tss_outdir (str): Path to time series data product output directory
        targ (str): Target name
    """

    for grat in ["g160m", "g230l"]:
        if not os.path.exists(tss_outdir):
            os.makedirs(tss_outdir)
        # First create the exposure level time-series spectra
        indir = os.path.join(datadir, grat, "exp")
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{VERSION.lower()}_tss.fits")
        timeseries.process_files(grat.upper(), outfile, indir, overwrite=True) 
        
        indir = os.path.join(datadir, grat, "split")
        outfile = os.path.join(tss_outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{VERSION.lower()}_split-tss.fits")
        timeseries.process_files(grat.upper(), outfile, indir,
                                 BINS[targ][grat]["wave"], BINS[targ][grat]["min_exptime"], overwrite=True)

    return tss_outdir


def main(datadir, orig_datadir, tss_outdir, targ, custom_caldir=None,
         g230l_disptab=G230L_DISPTAB):
    """
    Perfroms the timeseries data product calibration and creation.

    Args:
        datadir (str):
        orig_datadir (str):
        tss_outdir (str):
        targ (str):
        custom_caldir (str):
        g230l_disptab (str):
    """

    copy_origdata(datadir, orig_datadir)
    calibrate_data(datadir, custom_caldir, g230l_disptab)
    copy_caldata(datadir)
    create_splittags(datadir, targ)
    correct_vignetting(datadir)
    create_timeseries(datadir, tss_outdir, targ)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="datadir",
                        help="Path to directory with timeseries data")
    parser.add_argument(dest="orig_datadir",
                        help="")
    parser.add_argument(dest="tss_outdir",
                        help="")
    parser.add_argument(targ="targ",
                        help="Name of target")
    parser.add_argument(dest="custom_caldir",
                        help="")
    parser.add_argument(dest="disptab", default=G230L_DISPTAB,
                        help="If specified, custom COS NUV DISPTAB to use for calibration")
    args = parser.parse_args()

    main(args.datadir, args.orig_datadir, args.tss_outdir, args.targ, args.custom_caldir, args.disptab)
