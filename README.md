# ULLYSES

This repo contains the codes used to create the high level science products (HLSPs) for the targets in the ULLYSES special program. See more info about ULLYSES and the targets here: ullyses.stsci.edu.

### Installation

JWST pipeline readme - has info on how to do conda
add link to conda/anaconda

clone the repo

python setup.py install

should also install ullyses-utils repo (link here)

### Creating HLSPs

Individual HLSPs for a single target can be created by putting all of the input
files into one directory.  The input files are _x1d.fits files for COS and STIS,
and _vo.fits for FUSE.  You can create all the HLSPs for this target by running
the wrapper script, which can be done from the command line.  For convenience,
it is recommended to create an environment variable pointing to the location
of the wrapper script:

    export ubin=/path/to/github/checkout/ullyses

Then invoke the script from the directory containing the files to be processed:

    cd /directory/containing/data/files/
    python $ubin/wrapper.py -o './products'

Alternatively, the script can be run from a directory that doesn't contain the data
to be processed by using the ``-i /directory/containing/data/`` option:

    python $ubin/wrapper.py -o /directory/to/put/products -i /directory/containing/input/data

splittag_wrapper.py


timeseries.py


### Contributing

If you want to suggest changes to this content do the following:

1. Fork it.
2. Create your feature branch (git checkout -b my-new-feature).
3. Add your changes to staging area (git add myfile); This can be repeated multiple times.
4. If you are adding a new style guide, do not forget to update guides listing at README.md.
5. Commit your changes in staging area (git commit -m 'Added some feature').
6. Push to the branch (git push origin my-new-feature).
7. Create new Pull Request (PR).
8. Ask for a PR review.
