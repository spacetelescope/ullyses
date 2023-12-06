import argparse
import datetime
import shutil
import numpy as np
import os
import glob
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as pl
from stistools import x1d
import subprocess

from ullyses.calibrate_stis_data import calibrate_stis_data
from ullyses.stis_coadd_x1d import coadd_1d_spectra
import ullyses_utils

utils_dir = ullyses_utils.__path__[0]
from ullyses_utils.readwrite_yaml import read_config
from . import __release__

CONFIG_DIR = os.path.join(utils_dir, "data", "stis_configs")
OUTDIR_ROOT = None
nowdt = datetime.datetime.now()
if OUTDIR_ROOT is None:
    OUTDIR_ROOT = nowdt.strftime("%Y%m%d_%H%M")

"""
This code is used to preform custom calibrations
 of STIS MAMA and CCD data.

Arguments:
    datadir (str): Path to working data directory
    targ (str): Target Name
    orig_datadir (str): Path to original raw data
    outdir (str): Directory for output calibrated data
    config_dir (str): Path to directory with STIS configuration files
    copydir (str): Name of directory to copy output products to

Examples of the YAML configuration file for STIS targets
are located in the ullyses-utils github repo:
https://github.com/spacetelescope/ullyses-utils/tree/main/utils/data/stis_configs.
"""


def copy_origfiles(datadir, orig_datadir):
    """
    Copy the original STIS files to ensure we leave them intact.

    Args:
        datadir (str): Path to directory to copy data to
        orig_datadir (str): Path to directory with original data

    Return:
        None
    """

    files = glob.glob(os.path.join(orig_datadir, "o*fits"))
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in files:
        shutil.copy(item, datadir)

    print(f"\nCopied TTS data from {orig_datadir} to {datadir}\n")


def copy_cvso104_fringeflat(cvso109_dir, cvso104_fringeflat):
    """
    A special fringeflat is needed for CVSO-109.
    Copy CVSO-104's fringeflat to use with CVSO-109.
    """

    if not os.path.exists(cvso109_dir):
        os.makedirs(cvso109_dir)
    shutil.copy(cvso104_fringeflat, cvso109_dir)


def make_custom_x1ds(datadir, outdir, targ, config_dir=CONFIG_DIR):
    """
    Calibrate the STIS x1ds using the calibrate_stis_data function.

    Args:
        datadir (str):
        outdir (str):
        targ (str):
        config_dir (str):

    Return:
        None
    """

    configs = glob.glob(os.path.join(config_dir, f"{targ}_*yaml"))
    assert len(configs) != 0, f"No config files for target {targ}"
    for config in configs:
        calibrate_stis_data(datadir, config, outdir=outdir)
    print(f"\nMade custom products for target {targ}, wrote to {outdir}\n")


def coadd_blended_spectra(x1ds, targ, outdir):
    """
    Coadd STIS x1ds that have blended spectra.

    VERY IMPORTANT: the second listed target for each set of two needs to be the
    companion. The DQ arrays of the companion are ignored with coadding the two
    spectra.

    Args:
        x1ds (list): List of x1ds to coadd
        targ (str): Target name
        outdir (str): Path to output directory

    Return:
        None
    """

    coadd_1d_spectra(x1ds, targ, outdir)
    print(f"\nMade coadded blended spectra for {targ}, wrote to {outdir}\n")


def copy_rename_yaml(targ, outdir, config_dir=CONFIG_DIR):
    """
    Copy the YAML file from the configuration directory
    to the output directory.

    Args:
        targ (str): Target name
        outdir (str): Path to output directory
        config_dir (str): Path to directory with configuration file

    Return:
        None
    """

    yamlfiles = glob.glob(os.path.join(config_dir, f"{targ}_*yaml"))
    for item in yamlfiles:
        f = os.path.basename(item)
        spl = f.split("_")
        targ = spl[0].lower()
        grating = spl[1].split(".")[0]
        newname = f"hlsp_ullyses_hst_stis_{targ}_{grating}_{__release__}_spec.yaml"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        shutil.copyfile(item, os.path.join(outdir, newname))
    print(f"\nCopied and renamed {targ} YAML files from {config_dir} to {outdir}\n")


def copy_products(files, copydir):
    """
    Once you have created custom STIS 1D spectra, you may
    copy them to another destination.

    Args:
        outdir (str): Path to output directory
        copydir (str): Path to copy output to

    Return:
        None
    """

    for item in files:
        if not os.path.exists(copydir):
            os.makedirs(copydir)
        shutil.copy(item, copydir)
    print(f"\nCopied custom x1ds from {outdir} to {copydir}\n")


def main(datadir, targ, orig_datadir, outdir, config_dir=CONFIG_DIR, copydir=None):
    """
    Run script to make custom x1d files for STIS.

    Args:
        datadir (str): Path to working data directory
        targ (str): Target Name
        orig_datadir (str): Path to original raw data
        outdir (str): Directory for output calibrated data
        config_dir (str): Path to directory with STIS configuration files
        copydir (str): Name of directory to copy output products to

    Return:
        None
    """

    copy_origfiles(datadir, orig_datadir, targ)
    # copy_cvso104_fringeflat(cvso109_dir, cvso104_fringeflat)
    make_custom_x1ds(datadir, outdir, targ, config_dir)
    # coadd_blended_spectra(x1ds, targ, outdir)
    copy_rename_yaml(targ, outdir, config_dir)
    if copydir is not None:
        copy_products(outdir, copydir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="datadir",
                        help="Path to working data directory")
    parser.add_argument(dest="orig_datadir", default=None,
                        help="Path to original raw data. This is copied to datadir to protect the originals")
    parser.add_argument(dest="outdir",
                        help="Directory for output calibrated data")
    parser.add_argument(dest="targ",
                        help="Name of target to calibrate")
    parser.add_argument("--config_dir", default=CONFIG_DIR,
                        help="Path to directory with STIS configuration files")
    parser.add_argument("--copydir", default=None,
                        help="Name of directory to copy output products to")

    args = parser.parse_args()
    main(args.datadir, args.targ, args.orig_datadir, args.outdir, args.config_dir, args.copydir)
