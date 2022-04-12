# ULLYSES

This repo contains the codes used to create the high level science products (HLSPs) for the targets in the ULLYSES special program. See more info about ULLYSES and the targets here: ullyses.stsci.edu.

### Installation

JWST pipeline readme - has info on how to do conda
add link to conda/anaconda

clone the repo

python setup.py install

should also install ullyses-utils repo (link here)

### Creating HLSPs

wrapper.py


splittag_wrapper.py 


timeseries.py
#### Creating custom-calibrated individual spectra
Prior to turning spectra into ULLYSES HLSPs, some targets require extra reprocessing to
fix various calibration issues. For example, STIS/G750L data need to be defringed, and
COS/G230L data need to be corrected to remove vignetting effects. Once these custom
calibration steps have been applied, a keyword is added to the output FITS file signifying
that the file should be considered a *level0* HLSP, or custom-calibrated *individual*
1D spectrum. The various level0 products, and how to create them, are described below.

**Custom-calibrated STIS spectra**
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
where products and diagnostic plots will be written.

Other optional arguments are:
```
  -c, --clobber                     If True, overwrite existing products
  --nolog                           If True, do not produce log file
  -l LOGFILE, --logfile LOGFILE
                                    Name of output log file
```

**TODO:** add info about coadding blended spectra.

**Wavelength-shifted COS spectra**
If a target is not perfectly centered in the COS aperture, the wavelength array can be
offset from its true values. Wavelength offsets can be easily corrected by recalibrating
the data and supplying a shift file to CalCOS, as described in the 
[COS Data Handbook](https://hst-docs.stsci.edu/cosdhb/chapter-5-cos-data-analysis/5-3-working-with-extracted-spectra#id-5.3WorkingwithExtractedSpectra-5.3.2RedoingSpectralExtraction).

The ULLYSES team has identified stars which require such wavelength offset corrections, and 
documented the necessary shifts in text files stored in the 
[ullyses-utils](https://github.com/spacetelescope/ullyses-utils/tree/main/utils/data/cos_shifts) 
repository.

You can easily create wavelength-shifted COS spectra using the pre-defined ULLYSES shift 
files like so:

Once you have a shift file, you can create wavelength-shifted COS spectra by operating on
a directory with multiple exposures, or by operating on a single file. You can easily 
create wavelength-shifted COS spectra using the pre-defined ULLYSES shift files like so
```
python apply_cos_shifts.py <infiledir> <outdir> 
```
where `<infiledir>` is the input filename or directory of files that should be shifted,
and `<outdir>` is the directory to write shifted 1D spectra.
In this case, the target name in the input file(s) header(s) *must* match the target name
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

**Custom-flagged FUSE spectra**
The ULLYSES team uses FUSE 
[VO (Virtual Observatory)](https://ui.adsabs.harvard.edu/abs/2007PASP..119..527D/abstract) 
files with minimal modification. A DQ (Data Quality) array is added to each VO file, as is
required by the ULLYSES pipeline. For the majority of FUSE targets, this DQ array is 
uniformly zero, meaning there are no data quality issues. However, for a handful
of targets, custom flagging is imposed in order to screen out bad spectral regions. 
Possible DQ flags include: 
* DQ=1 (region was affected by the worm) 
* DQ=2 (poor photometric quality) 

To create these custom-flagged FUSE spectra:
**Still need to fill this in**

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
