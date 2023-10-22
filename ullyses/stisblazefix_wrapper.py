"""
Apply the stisblazefix to selected echelle data.
"""

from stisblazefix import fluxfix
import pandas as pd
from astropy.io import fits
import shutil
import ullyses_utils


def apply_sbf(infile, outdir=None, pdfname='/dev/null', verbose=True):
    fluxfix([infile], pdfname=pdfname)
    newfile = infile.replace("x1d.fits", "x1f.fits")
    if outdir is not None:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        newfile = os.path.join(outdir, os.path.basename(newfile))
        shutil.move(newfile, outfile)
    else:
        outfile = newfile
    with fits.open(outfile, mode="update") as hdulist:
    # Since custom processing was performed, mark these as level0
        hdulist[0].header["HLSP_LVL"] = 0
    if verbose is True:
        print(f"Wrote {outfile}")


