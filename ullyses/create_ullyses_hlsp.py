import re
import glob
import os
import astropy
import pandas as pd
from astropy.io import fits
import numpy as np
import datetime
from datetime import datetime as dt
from astropy.time import Time

from ullyses.ullyses_hlsp import Ullyses
import ullyses_utils
from ullyses_utils.ullyses_config import VERSION, CAL_VER, RENAME
from ullyses_utils.parse_csv import parse_aliases

UTILS_DIR = ullyses_utils.__path__[0]
PHOT_DIR = os.path.join(UTILS_DIR, "data", "lcogt_photometry")
RED = "\033[1;31m"
RESET = "\033[0;0m"


def make_imaging_hlsps(drcfile, outdir, targ, hdr_targ=None, hlspname=None):
    targ = targ.lower()
    if hdr_targ is None:
        hdr_targ = targ
    filt = fits.getval(drcfile, "filter")
    ra = fits.getval(drcfile, "ra_targ")
    dec = fits.getval(drcfile, "dec_targ")
    if hlspname is None:
        hlspname = f"hlsp_ullyses_hst_wfc3_{targ}_{filt.lower()}_{VERSION}_drc.fits"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, hlspname)
    U = Ullyses(files=[drcfile], hlspname=outfile, targname=hdr_targ, ra=ra, dec=dec, 
                cal_ver=CAL_VER, version=VERSION, level=6, hlsp_type="drizzled")
    U.make_hdrs_and_prov()
    U.write_file()


def make_lcogt_tss(indir, outdir, targ, hlspname=None, photfile=None):
    if photfile is None:
        photfile = os.path.join(PHOT_DIR, f"{targ.upper()}_phot.txt")
    aliases = parse_aliases()
    alias_mask = aliases.apply(lambda row: row.astype(str).str.fullmatch(re.escape(targ.upper())).any(), axis=1)
    if set(alias_mask) != {False}:
        ull_targname = aliases[alias_mask]["ULL_MAST_name"].values[0]
        targetinfo_file = os.path.join(UTILS_DIR, "data", "target_metadata", "pd_targetinfo.json")
        master_list = pd.read_json(targetinfo_file, orient="split")
        master_list = master_list.apply(lambda x: x.astype(str).str.upper())
        try:
            coords = master_list.loc[master_list["mast_targname"] == ull_targname][["ra", "dec"]].values
            ra, dec = coords[0]
        except:
            print(f"{RED}NO COORDINATES FOUND FOR {ull_targname}{RESET}")
            ra, dec = (0, 0)
    else:
        ull_targname = targ
        ra, dec = (0, 0)
        print(f"{RED}NO ALIAS AND COORDINATES FOUND FOR {ull_targname}{RESET}")

    df = pd.read_csv(photfile, 
            names=["filename", "mjdstart", "mjdend", "wl", "flux", "err"],
            skiprows=[0], delim_whitespace=True)
    df = df.sort_values("mjdstart")
    filenames = df.filename.tolist()
    filepaths = [os.path.join(indir, x) for x in filenames]
    if hlspname is None:
        if 3560 in df["wl"].values:
            filts = "uprime-v-iprime"
        else:
            filts = "v-iprime"
        file_targname = rename_target(ull_targname)
        hlspname = f"hlsp_ullyses_lcogt_04m_{file_targname.lower()}_{filts}_{VERSION}_tss.fits"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, hlspname)
    U = Ullyses(files=filepaths, hlspname=outfile, targname=ull_targname, 
                ra=ra, dec=dec, cal_ver=CAL_VER, version=VERSION, level=5, 
                hlsp_type="lcogt", photfile=photfile)
    U.make_hdrs_and_prov()
    U.write_file()


def rename_target(targname):
    targ_lower = targname.lower()
    if targ_lower in RENAME:
        file_targname = RENAME[targ_lower]
    else:
        file_targname = targname
    return file_targname
