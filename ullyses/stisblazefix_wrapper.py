"""
Apply the stisblazefix to selected echelle data.
"""

from stisblazefix import fluxfix
import pandas as pd

import ullyses_utils


def apply_sbf(filename, pdfname = '/dev/null'):
    fluxfix([filename], pdfname=pdfname)


def ullyses_sbf():
    ull_dbfile = ullyses_utils.__path__[0] + '/data/target_metadata/ullyses_db.csv'
    db = pd.read_csv(ull_db)
    files = db.loc[db["sbf"] == True]
    for item in files:
        apply_sbf(item)

