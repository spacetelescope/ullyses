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

from readwrite_yaml import read_config
from calibrate_stis_data import wrapper
from stis_coadd_x1d import coadd_1d_spectra

TARGS = [# DR2
       "cvso-104", "cvso-107", "cvso-109", "cvso-146", "cvso-165", "cvso-17",
       "cvso-176", "cvso-36", "cvso-58", "cvso-90", "v-tx-ori", "v505-ori",
       "v510-ori", 
       #dr3
       "sz10", "sz45", "sz69", "sz71", "sz72", "sz75", "sz77", "v-in-cha",
       "v-xx-cha", "chx18n", "2massj11432669-7804454", "echa-j0844.2-7833",
       "hn5", "sstc2dj160000.6-422158"]

VERSION = "dr3"
HLSPDIR = "/astro/ullyses/ULLYSES_HLSP"
#HLSPDIR = "ULLYSES_HLSP" 
VETTED_DIR = "/astro/ullyses/all_vetted_data_dr3"
#VETTED_DIR = "all_vetted_data_dr3"
DATADIR = "/astro/ullyses/ULLYSES_DATA"
CUSTOM_CAL = "/astro/ullyses/custom_cal"
CONFIG_DIR = "config_files/"
#CONFIG_DIR = "test_configs/"
OUTDIR_ROOT = None
nowdt = datetime.datetime.now()
if OUTDIR_ROOT is None:
    OUTDIR_ROOT = nowdt.strftime("%Y%m%d_%H%M")

def copy_rawfiles():
    for targ in TARGS:
        files = glob.glob(os.path.join(DATADIR, targ.upper(), "o*fits"))
        # Targs with periods in their name must be specially renamed or defringe will crash
        if "." in targ:
            spl = targ.split(".")
            targ = "".join(spl)
        destdir = os.path.join(CUSTOM_CAL, targ.lower())
        if not os.path.exists(destdir):
            os.makedirs(destdir)
        for item in files:
            shutil.copy(item, destdir)

    # A special fringeflat is needed for CVSO-109, copy it
    origin = os.path.join(DATADIR, "CVSO-104", "oe9k1s050_raw.fits")
    dest = os.path.join(CUSTOM_CAL, "cvso-109")
    shutil.copy(origin, dest)

    print(f"\nCopied TTS data from {DATADIR} to {CUSTOM_CAL}\n")


def copy_products(outdir_root=OUTDIR_ROOT):
    for targ in TARGS:
        if "." in targ:
            spl = targ.split(".")
            dirtarg = "".join(spl)
        else:
            dirtarg = targ
        outdir = os.path.join(CUSTOM_CAL, dirtarg.lower(), outdir_root)
        files = glob.glob(os.path.join(outdir, "*x1d.fits"))
        for item in files:
            actualtarg = fits.getval(item, "targname").lower()
            destdir = os.path.join(VETTED_DIR, actualtarg)
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            shutil.copy(item, destdir)
    print(f"\nCopied TTS final products from {CUSTOM_CAL}/*/{outdir_root} to {VETTED_DIR}\n")


def make_custom_x1ds(outdir_root=OUTDIR_ROOT):
    configs = glob.glob(os.path.join(CONFIG_DIR, "*yaml"))
    for config in configs:
        config_name = os.path.basename(config)
        if config_name == "target_grating.yaml":
            continue
        targ = config_name.split("_")[0].lower()
        if "." in targ:
            spl = targ.split(".")
            dirtarg = "".join(spl)
        else:
            dirtarg = targ
        indir = os.path.join(CUSTOM_CAL, dirtarg)
        outdir = os.path.join(indir, outdir_root)
        wrapper(indir, config, outdir=outdir)
    print(f"\nMade custom products for data in {CUSTOM_CAL}, wrote to {CUSTOM_CAL}/*/{outdir_root}\n")


def coadd_blended_spectra(outdir_root=OUTDIR_ROOT):
    d = {"cvso-109": [["oe9k2s020_cvso-109a_x1d.fits", "oe9k2s020_cvso-109b_x1d.fits"],
                      ["oe9k2s030_cvso-109a_x1d.fits", "oe9k2s030_cvso-109b_x1d.fits"],
                      ["oe9k2s010_cvso-109a_x1d.fits", "oe9k2s010_cvso-109b_x1d.fits"]],
         "cvso-165": [["oe9j2s010_cvso-165a_x1d.fits", "oe9j2s010_cvso-165b_x1d.fits"]]}
    for targ in d:
        for pair in d[targ]:
            files0 = pair
            indir = os.path.join(CUSTOM_CAL, targ.lower(), outdir_root)
            files = [os.path.join(indir, x) for x in files0]
            coadd_1d_spectra(files, targ, outdir=indir)
    print(f"\nMade coadded blended spectra for {d.keys()}, wrote to {CUSTOM_CAL}/*{outdir_root}\n")

def copy_rename_yaml():
    yamlfiles0 = glob.glob(os.path.join(CONFIG_DIR, "*yaml"))
    yamlfiles = [x for x in yamlfiles0 if os.path.basename(x) != "target_grating.yaml"]
    for item in yamlfiles:
        f = os.path.basename(item)
        spl = f.split("_")
        targ = spl[0].lower()
        grating = spl[1].split(".")[0]
        newname = f"hlsp_ullyses_hst_stis_{targ}_{grating}_{VERSION}_spec.yaml"
        dest = os.path.join(HLSPDIR, targ, VERSION)
        if not os.path.exists(dest):
            os.makedirs(dest)
        shutil.copyfile(item, os.path.join(dest, newname))
    print(f"\nCopied and renamed YAML files from {CONFIG_DIR} to {HLSPDIR}\n")
    

def main():
#    copy_rawfiles()
    make_custom_x1ds()
    coadd_blended_spectra()
    copy_rename_yaml()
    copy_products()


if __name__ == "__main__":
    main()
