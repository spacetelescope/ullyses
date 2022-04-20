# ULLYSES

This repo contains the codes used to create the high level science products (HLSPs) for the targets in the Hubble Space Telescope’s (HST) Ultraviolet Legacy Library of Young Stars as Essential Standards (ULLYSES) program. See more info about ULLYSES and its targets at [ullyses.stsci.edu](https://ullyses.stsci.edu).

A full description of the data products produced by the ULLYSES team can be found at [ULLYSES Data Products](https://ullyses.stsci.edu/ullyses-data-description.html). 

### Installation

The `ullyses` package can be installed into a virtualenv or conda environment via `pip`. We recommend that for each installation you start by creating a fresh environment that only has Python installed and then install the `ullyses` package and its dependencies into that bare environment. If using conda environments, first make sure you have a recent version of Anaconda or Miniconda installed.

The first two steps are to create and activate an environment:

    conda create -n <env_name> python=3.8
    conda activate <env_name>
   
Python version 3.8 or greater is required for some dependencies, including `calcos`, the COS data calibration pipeline used in these scripts.

To install your own copy of the code into that environment, you first need to fork and clone the `ullyses` repo:

    cd <where you want to put the repo>
    git clone https://github.com/spacetelescope/ullyses
    cd ullyses
    
*Note: `python setup.py install` and `python setup.py develop` commands do not work.*

Install from your local checked-out copy:

    pip install .

If you anticipate making edits to any of the ULLYSES scripts, it may be more useful to 
install from your local checked-out copy as an "editable" install:

    pip install -e .

All package dependencies will be installed simultaneously, including `ullyses-utils`, which can be found at https://github.com/spacetelescope/ullyses-utils.

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


### Contributions and Feedback

We welcome contributions and feedback on this project. If you want to suggest changes to this content, please do the following:

1. Fork it.
2. Create your feature branch (git checkout -b my-new-feature).
3. Add your changes to staging area (git add myfile); This can be repeated multiple times.
4. If you are adding a new style guide, do not forget to update guides listing at README.md.
5. Commit your changes in staging area (git commit -m 'Added some feature').
6. Push to the branch (git push origin my-new-feature).
7. Create new Pull Request (PR).
8. Ask for a PR review.

We strive to provide a welcoming community to all of our users by abiding with
the [Code of Conduct](CODE_OF_CONDUCT.md).

If you have questions or concerns regarding the software, please open an issue at https://github.com/spacetelescope/ullyses/issues or contact the [HST Help Desk](https://hsthelp.stsci.edu). If you have questions regarding the ULLYSES program design or data, please contact the [HST Help Desk](https://hsthelp.stsci.edu).

