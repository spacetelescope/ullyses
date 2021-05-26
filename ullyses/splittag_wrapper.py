import os
import calcos
from costools import splittag
import glob

"""
This code will:
(1) use costools.splittag to divide corrtag files for the variable stars into multiple corrtags based on a time interval
(2) run calcos on the output corrtags to re-extract them into x1ds (x1dsums?)
(3) remove unnecessary file outputs such as counts images to maintain a tidy directory

"""

DATADIR = '/astro/ullyses/ULLYSES_DATA/'
OUTPUTDIR = '/astro/ullyses/'




splittag.splittag("rootname_corrtag_a.fits", "split250s",
                  starttime=None, increment=None, endtime=None,
                  time_list="0, 250, 500, 750, 1000")



if __name__ == "__main__":

    # need to get the names of all the targets we will be making lightcurves for
    # probably a better way to do this than hard-coding. will revisit
    targetnamelist = ['TW-HYA']

    for target in targetnamelist:

        targetdir = os.path.join(DATADIR, target)

        corrtags = glob.glob(os.path.join(targetdir, '*corrtag_a*.fits'))

