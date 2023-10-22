"""
Apply the stisblazefix to selected echelle data.
"""

from stisblazefix import fluxfix
import pandas as pd
from astropy.io import fits
import shutil
import ullyses_utils


def apply_sbf(infile, outfile=None, pdfname='/dev/null', verbose=True):
    fluxfix([infile], pdfname=pdfname)
    newfile = infile.replace("x1d.fits", "x1f.fits")
    if outfile is not None:
        shutil.move(newfile, outfile)
    else:
        outfile = newfile
    with fits.open(outfile, mode="update") as hdulist:
    # Since custom processing was performed, mark these as level0
        hdulist[0].header["HLSP_LVL"] = 0
    if verbose is True:
        print(f"Wrote {outfile}")


