"""
Apply the stisblazefix to selected echelle data.
"""

from stisblazefix import fluxfix
import pandas as pd

import ullyses_utils


def apply_sbf(filename, pdfname = '/dev/null', verbose=True):
    fluxfix([filename], pdfname=pdfname)
    if verbose is True:
        newfile = filename.replace("x1d.fits", "x1f.fits")
        print(f"Wrote {newfile}")


def ullyses_sbf():
    print("Running stisblazefix...")
    ull_dbfile = ullyses_utils.__path__[0] + '/data/target_metadata/ullyses_db.csv'
    db = pd.read_csv(ull_db)
    files = db.loc[db["sbf"] == True]
    for item in files:
        apply_sbf(item)

