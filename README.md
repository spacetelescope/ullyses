# ULLYSES

This repo contains the codes used to create the high level science products (HLSPs) for the targets in the Hubble Space Telescope’s (HST) Ultraviolet Legacy Library of Young Stars as Essential Standards (ULLYSES) program. See more info about ULLYSES and its targets at [ullyses.stsci.edu](https://ullyses.stsci.edu).

A full description of the data products produced by the ULLYSES team can be found at [ULLYSES Data Products](https://ullyses.stsci.edu/ullyses-data-description.html). 

## Installation

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

## Creating HLSPs

There are four main types of HLSPs that ULLYSES produces:
1. coadded spectra
2. timeseries spectra
3. re-packaged drizzled images
4. custom-calibrated individual spectra, called "level0" HLSPs

Below are instructions for creating each type of HLSP.

### Coadded Spectra HLSPs 
Individual coadded products for a single target can be created by putting all of the input
files into one directory. The input files are `_x1d.fits` files for COS and STIS,
and `_vo.fits` for FUSE. These input files may also be level0 (custom-calibrated spectra) 
themselves. You can create all the coadd HLSPs for this target by running
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

### Timeseries Spectra HLSPs

There are two main flavors of timeseries spectra: exposure level and sub-exposure level.
Exposure level timeseries spectra are essentially stacked individual 1D spectra.
Sub-exposure level timeseries spectra are made from _split_ 1D spectra. That is, for 
time-tag data (COS and STIS UV data), exposures are broken down into even smaller time 
chunks. 

#### HST Timeseries

For HST, both exposure and sub-exposure timeseries spectra are created the same way.
To create a timeseries spectrum, you must supply a configuration YAML file. The ULLYSES
team has already created such files for all monitoring stars and recorded the optimal
parameters in YAML files stored in the
[ullyses-utils](https://github.com/spacetelescope/ullyses-utils/tree/main/utils/data/timeseries)
repository. You may use the ULLYSES YAML files as input, or supply your own, but they
must conform to the format outlined here. TODO- ADD LINK TO TEMPLATE FILE

Once you have a YAML file, you create the timeseries spectrum like so:
```
python ctts_cal.py --orig <origdir> --copydir <copydir> --hlspdir <hlspdir> --targ <targ> --yaml <yaml>
```

where `<origdir>` is the directory which houses all the input data, `<copydir>` is the directory
to copy input data to (data will be modified), `<hlspdir>` is the directory to write the final
timeseries spectra, `<targ>` is the ULLYSES name of the target being calibrated, and `<yaml>` is the
YAML confirmation file.

#### Photometric Timeseries

It is possible to create a timeseries "spectrum" using photometric measurements over time.
The ULLYSES team has performed photometry on ULLYSES low-mass stars using the LCOGT network
of telescopes.

To create LCOGT photometric timeseries spectra:
```
python lcogt_hlsps_wrapper.py -i <indir> -o <outdir> -t <targ>
```

where `<indir>` is the directory which contain the original LCOGT FITS images (used to 
extract observational metadata), `<outdir>` is the directory to write the HLSPs, and
`<targ>` is the ULLYSES target name.

### Re-packaged Drizzled Image HLSPs

The ULLYSES team creates drizzled WFC3 images for the low-metallicity galaxies NGC3109 and
SextansA. These images are repackaged to conform to the ULLYSES HLSP requirements, but the 
data array values are left untouched. 

To create drizzled image HLSPs:
```
python imaging_hlsps_wrapper.py <drcfile> -o <outdir> -t <targ>
```

where `<drcfile>` is the name of the original drizzled DRC file, `<outdir>` is the 
directory to write the HLSP to, and `<targ>` is the ULLYSES target name.

Optional arguments are:
```
  --hdr_targ HDR_TARG   If specified, alternative target name to use in HLSP file name
  --hlspname HLSPNAME   Name of output HLSP file. By default, follows ULLYSES standard
```

### Custom-calibrated Spectra (level0) HLSPs

Prior to turning spectra into ULLYSES HLSPs, some targets require extra processing to
fix various calibration issues. For example, STIS/G750L data must be defringed or
wavelength offsets due to miscentering in the aperture must be corrected. Once these custom
calibration steps have been applied, a keyword is added to the output FITS file signifying
that the file should be considered a level0 HLSP- that is, a custom-calibrated *individual*
1D spectrum. The various level0 products, and how to create them, are described below.

#### Custom-calibrated STIS spectra
All T Tauri star STIS CCD observations, and a subset of STIS NUV- and FUV-MAMA observations,
require tailored calibrations. Special calibration steps can include: 
custom hot pixel identification and flagging, defringing for G750L observations, and 
customized spectral extraction parameters for T Tauri stars and their companions.

To create a custom-calibrated STIS spectra, you must supply a configuration YAML file.
The ULLYSES team has already examined each T Tauri star and recorded the optimal 
custom calibration parameters in YAML files stored in the 
[ullyses-utils](https://github.com/spacetelescope/ullyses-utils/tree/main/utils/data/stis_configs) 
repository. You may use the ULLYSES YAML files as input, or supply your own, but it must
conform to the 
[format outlined here](https://github.com/spacetelescope/ullyses-utils/blob/main/utils/data/stis_configs/target_grating.yaml).

Once you have a YAML file, you create the custom-calibrated STIS spectrum like so:
```
python calibrate_stis_data.py -i <indir> -y <yaml> -o <outdir>
```

where `<indir>` is the input directory that houses the data you wish to calibrate,
`<yaml>` is the name of the yaml file, and `<outdir>` is the output directory
where products, logs, and diagnostic plots will be written. A log file of the format
`YYYYMMDD_HHmm_cal.log` will be written, unless otherwise specified using the optional
arguments below.

Optional arguments are:
```
  -c, --clobber                     If True, overwrite existing products
  --nolog                           If True, do not produce log file
  -l LOGFILE, --logfile LOGFILE
                                    Alternative name of output log file
```

**TODO:** add info about coadding blended spectra.

#### Wavelength-shifted COS spectra
If a target is not perfectly centered in the COS aperture, the wavelength array can be
offset from its true values. Wavelength offsets can be easily corrected by recalibrating
the data and supplying a shift file to CalCOS, as described in the 
[COS Data Handbook](https://hst-docs.stsci.edu/cosdhb/chapter-5-cos-data-analysis/5-3-working-with-extracted-spectra#id-5.3WorkingwithExtractedSpectra-5.3.2RedoingSpectralExtraction).

The ULLYSES team has identified stars which require such wavelength offset corrections and 
documented the necessary shifts in text files stored in the 
[ullyses-utils](https://github.com/spacetelescope/ullyses-utils/tree/main/utils/data/cos_shifts) 
repository.

Once you have a shift file, you can create wavelength-shifted COS spectra by supplying on
a directory with multiple exposures, or by supplying on a single file. You can easily 
create wavelength-shifted COS spectra using the pre-defined ULLYSES shift files like so:
```
python apply_cos_shifts.py <infiledir> <outdir> 
```
where `<infiledir>` is the input filename or directory of files that should be shifted,
and `<outdir>` is the directory to write shifted 1D spectra.
In this case, the target name in the input file(s) header(s) **_must_** match the target name
in the shift file. If for some reason it does not, you must also supply the target name 
as it appears in the shift file using an additional `-t <targ>` argument.

You may also create wavelength-shifted COS spectra using your own custom shift file, which
is done like so:

```
python apply_cos_shifts.py <infiledir> <outdir> -s <shift_file>
```
where `<shift_file>` is the name of your custom shift file. 

Other optional arguments are:
```
  --copydir COPYDIR     Name of directory to copy shifted products to
  -c, --overwrite       If True, overwrite existing products
```

#### Custom-flagged FUSE spectra
The ULLYSES team typically uses FUSE 
[VO (Virtual Observatory)](https://ui.adsabs.harvard.edu/abs/2007PASP..119..527D/abstract) 
files with minimal modification. A DQ (Data Quality) array is added to each VO file, as is
required by the ULLYSES pipeline. For the majority of FUSE targets, this DQ array is 
uniformly zero, meaning there are no data quality issues. However, for a handful
of targets, custom flagging is imposed in order to screen out bad spectral regions. 
Possible DQ flags include: 
* DQ=1 (region was affected by the worm) 
* DQ=2 (poor photometric quality) 

To create these custom-flagged FUSE spectra:
```
python make_flagged_fuse.py -i <infile> -o <outdir>
```
where `<infile>` is the file to flag and `<outdir>` is the output directory.

Other optional arguments are:
```
  -t TARG or --targ TARG  ULLYSES name of target
  -c or --overwrite       If True, overwrite existing products
```

## Contributions and Feedback

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

