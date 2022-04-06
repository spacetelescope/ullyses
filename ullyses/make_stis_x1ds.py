import datetime
import shutil
import numpy as np
import os
import glob
import argparse

from calibrate_stis_data import calibrate_stis_data
from stis_coadd_x1d import coadd_1d_spectra
import ullyses_utils
utils_dir = ullyses_utils.__path__[0]
from ullyses_utils.ullyses_config import VERSION, RENAME
from ullyses_utils.readwrite_yaml import read_config

CONFIG_DIR = os.path.join(utils_dir, "data", "stis_configs")
OUTDIR_ROOT = None
nowdt = datetime.datetime.now()
if OUTDIR_ROOT is None:
    OUTDIR_ROOT = nowdt.strftime("%Y%m%d_%H%M")


def copy_origfiles(targ, datadir, orig_datadir):
    """
    Copy the original STIS files to ensure we leave them intact.
    """
    files = glob.glob(os.path.join(orig_datadir, "o*fits"))
    # Targs with periods in their name must be specially renamed or defringe will crash
    if "." in targ:
        assert targ in RENAME, f"Renaming scheme not known for {targ}"
        targ = RENAME[targ]
    if not os.path.exists(datadir):
        os.makedirs(datadir)
    for item in files:
        shutil.copy(item, datadir)

    print(f"\nCopied TTS data from {orig_datadir} to {datadir}\n")


def copy_cvso104_fringeflat(cvso109_dir, cvso104_fringeflat):
    """
    A special fringeflat is needed for CVSO-109.
    Copy CVSO-104's fringeflat to use with CVSO-109
    """
    if not os.path.exists(cvso109_dir):
        os.makedirs(cvso109_dir)
    shutil.copy(cvso104_fringeflat, cvso109_dir)


def make_custom_x1ds(datadir, outdir, targ, config_dir=CONFIG_DIR):
    """

    :param datadir:
    :param outdir:
    :param targ:
    :param config_dir:
    :return:
    """
    configs = glob.glob(os.path.join(config_dir, f"{targ}*yaml"))
    assert len(configs) != 0, f"No config files for target {targ}"
    if "." in targ:
        assert targ in RENAME, f"Renaming scheme not known for {targ}"
    calibrate_stis_data(datadir, configs, outdir)
    print(f"\nMade custom products for target {targ}, wrote to {outdir}\n")


def coadd_blended_spectra(x1ds, targ, outdir):
    """
    VERY IMPORTANT: the second listed target for each set of two needs to be the
    companion. The DQ arrays of the companion are ignored with coadding the two
    spectra. 
    """
    coadd_1d_spectra(x1ds, targ, outdir)
    print(f"\nMade coadded blended spectra for {targ}, wrote to {outdir}\n")


def copy_rename_yaml(targ, outdir, config_dir=CONFIG_DIR):
    yamlfiles = glob.glob(os.path.join(config_dir, f"{targ}*yaml"))
    for item in yamlfiles:
        f = os.path.basename(item)
        spl = f.split("_")
        targ = spl[0].lower()
        # Targs with periods in their name must be specially renamed
        if "." in targ:
            assert targ in RENAME, f"Renaming scheme not known for {targ}"
            targ = RENAME[targ]
        grating = spl[1].split(".")[0]
        newname = f"hlsp_ullyses_hst_stis_{targ}_{grating}_{VERSION}_spec.yaml"
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        shutil.copyfile(item, os.path.join(outdir, newname))
    print(f"\nCopied and renamed {targ} YAML files from {config_dir} to {outdir}\n")
    

def copy_products(outdir, copydir):
    """
    Once you have created custom STIS 1D spectra, you may
    copy them to another destination.
    """
    files = glob.glob(os.path.join(outdir, "*x1d.fits"))
    for item in files:
        if not os.path.exists(copydir):
            os.makedirs(copydir)
        shutil.copy(item, copydir)
    print(f"\nCopied custom x1ds from {outdir} to {copydir}\n")


def main(datadir, targ, orig_datadir, outdir, config_dir=CONFIG_DIR, copydir=None):
    copy_origfiles(datadir, orig_datadir, targ)
    #copy_cvso104_fringeflat(cvso109_dir, cvso104_fringeflat)
    make_custom_x1ds(datadir, outdir, targ, config_dir)
    #coadd_blended_spectra(x1ds, targ, outdir)
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
    parser.add_argument(dest= "targ",
                        help="Name of target to calibrate")
    parser.add_argument("--config_dir", default=CONFIG_DIR,
                        help="Path to directory with STIS configuration files")
    parser.add_argument("--copydir", default=None,
                        help="Name of directory to copy output products to")

    args = parser.parse_args()
    main(args.datadir, args.targ, args.orig_datadir, args.outdir, args.config_dir, args.copydir)
