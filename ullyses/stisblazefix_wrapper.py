"""
Apply the stisblazefix to selected echelle data.
"""

from stisblazefix import fluxfix
import pandas as pd
from astropy.io import fits

import ullyses_utils


def apply_sbf(filename, pdfname = '/dev/null', verbose=True):
    fluxfix([filename], pdfname=pdfname)
    newfile = filename.replace("x1d.fits", "x1f.fits")
    with fits.open(newfile, mode="update") as hdulist:
    # Since custom processing was performed, mark these as level0
        hdulist[0].header["HLSP_LVL"] = 0
    if verbose is True:
        print(f"Wrote {newfile}")


def ullyses_sbf():
    print("Running stisblazefix...")
    ull_dbfile = ullyses_utils.__path__[0] + '/data/target_metadata/ullyses_db.csv'
    db = pd.read_csv(ull_db)
    files = db.loc[db["sbf"] == True]
    for item in files:
        apply_sbf(item)

