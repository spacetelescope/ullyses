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
from ullyses_utils.parse_csv import parse_aliases
from ullyses_utils import match_aliases
from . import __version__, __release__

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
        hlspname = f"hlsp_ullyses_hst_wfc3_{targ}_{filt.lower()}_{__release__}_drc.fits"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, hlspname)
    U = Ullyses(files=[drcfile], hlspname=outfile, targname=hdr_targ, 
                cal_ver=__version__, version=__release__, level=6, hlsp_type="drizzled")
    U.make_hdrs_and_prov()
    U.write_file()


def make_lcogt_tss(indir, outdir, targ, hlspname=None, photfile=None):
    if photfile is None:
        photfile = os.path.join(PHOT_DIR, f"{targ.upper()}_phot.txt")
    hlsp_targname = match_aliases.match_aliases(targ)

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
        hlspname = f"hlsp_ullyses_lcogt_04m_{hlsp_targname.lower()}_{filts}_{__release__}_tss.fits"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, hlspname)
    ull_targname = match_aliases.match_aliases(hlsp_targname, "target_name_ullyses")
    U = Ullyses(files=filepaths, hlspname=outfile, targname=ull_targname, 
                cal_ver=__version__, version=__release__, level=5,
                hlsp_type="lcogt", photfile=photfile)
    U.make_hdrs_and_prov()
    U.write_file()

def make_xsu_hlsps(infile, outdir, targ, hlspname=None):
    hlsp_targname = match_aliases.match_aliases(targ)

    if hlspname is None:
        hlspname = f"hlsp_ullyses_vlt_xshooter_{hlsp_targname.lower()}_uvb-vis_{__release__}_vltspec.fits"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = os.path.join(outdir, hlspname)
    ull_targname = match_aliases.match_aliases(hlsp_targname, "target_name_ullyses")
    U = Ullyses(files=[infile], hlspname=outfile, targname=ull_targname,  
                cal_ver=__version__, version=__release__, level=7,
                hlsp_type="xsu", overwrite=True)
    U.make_hdrs_and_prov()
    U.write_file()
