"""
Apply the stisblazefix to selected echelle data.
"""

import os
import sys
try:
    from stisblazefix import fluxfix
except ModuleNotFoundError:
    print("To download stisblazefix, follow instructions here: https://stisblazefix.readthedocs.io/en/latest/#installation")
    sys.exit()
import pandas as pd
from astropy.io import fits
import shutil
import ullyses_utils


def apply_sbf(infile, outfile=None, pdfname='/dev/null', verbose=True):
    fluxfix([infile], pdfname=pdfname)
    x1ffile = infile.replace("x1d.fits", "x1f.fits")
    if outfile is None:
        outfile = x1ffile
    else:
        outdir = os.path.dirname(outfile)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        shutil.move(x1ffile, outfile)
    # Since custom processing was performed, mark these as level0
    with fits.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    if verbose is True:
        print(f"Wrote {outfile}")
    return outfile


