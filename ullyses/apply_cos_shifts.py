import argparse
from astropy.io import fits
import calcos
import os
import shutil
import glob
import pandas as pd
import datetime
import ullyses_utils
utils_dir = ullyses_utils.__path__[0]


"""
This code is used to apply custom wavelength shifts 
for ULLYSES targets that have wavelength calibration issues.

Users will need to create a text file with the columns
"rootname", "fppos", "flash", "segment", and "shift1"
for each target that needs correction.

rootname: the unique rootname for each rawtag to be corrected
FPPOS: typical input value is "any"
flash: typical input value is "any"
segment: FUVA, FUVB, NUVA, NUVB, or NUVC
shift1: the shift in x pixels to apply for that rawtag

Users will need to specify the location of the data 

"""

def apply_shifts_file(infile, outdir, shift_file):

    """
    Apply shifts from input txt file. The txt file is
    given as an input to calcos, which recalibrates the
    rawtag files.
    :return:
    """
    
    infile_name = os.path.basename(infile)
    if not infile_name.endswith("_asn.fits"):
        raise TypeError("Input file must be an association file")

    df = pd.read_csv(shift_file, delim_whitespace=True,
                     names=["rootname", "fppos", "flash", "segment", "shift1"],
                     header=None)

    ipppss = os.path.basename(infile)[:6]
    if ipppss not in df["rootname"].values:
        raise RuntimeError(f"Specified input file {infile} not in shifts file")
        
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    calcos.calcos(infile, shift_file=shift_file, outdir=outdir, 
                  verbosity=0)


def apply_shifts_dir(indir, outdir, shift_file):

    """
    Apply shifts from input txt file. The txt file is
    given as an input to calcos, which recalibrates the
    rawtag files.
    :return:
    """

    seen = []                                                                                                   
    df = pd.read_csv(shift_file, delim_whitespace=True,                                                         
                     names=["rootname", "fppos", "flash", "segment", "shift1"],                                 
                     header=None)
    uniq_ipppss = list(set(df["rootname"]))
    for ipppss in uniq_ipppss:
        files = glob.glob(os.path.join(indir, f"{ipppss}*asn.fits"))
        if  len(files) == 0:
            raise FileNotFoundError(f"{ipppss}* files are in shift file, but not in specified directory")

        apply_shifts_file(files[0], outdir, shift_file)


def determine_file_shifts(infile, targ=None):
    if targ is None:
        targ = fits.getval(infile, "targname")
    targ = targ.lower()
    shift_file = os.path.join(utils_dir, f"data/cos_shifts/{targ}_shifts.txt")
    if not os.path.exists(shift_file):
        raise FileNotFoundError(f"Shift file {shift_file} not found")

    return shift_file

def determine_dir_shifts(indir, targ=None):
    all_files = glob.glob(os.path.join(indir, targ, "*asn.fits"))
    if targ is not None:
        targ = targ.lower()
        shift_file = determine_file_shifts(all_files[0], targ=targ)
    else:
        shift_files = []
        for item in all_files:
            shift_files.append(determine_file_shifts(item))
        if len(list(set(shift_files))) != 1:
            raise ValueError(f"More than one target found in input directory {indir}")
        shift_file = shift_files[0]

    return shift_file


def add_hlsp_lvl0(outdir):
    x1ds = glob.glob(os.path.join(outdir, "*x1d.fits"))
    for x1d in x1ds:
        with fits.open(x1d, mode="update") as hdulist:
            hdulist[0].header["HLSP_LVL"] = 0
    

def copy_output_x1ds(outdir, copydir, overwrite=False):
    x1ds = glob.glob(os.path.join(outdir, "*x1d.fits"))
    if not os.path.exists(copydir):
        os.makedirs(copydir)

    for x1d in x1ds:
        if overwrite is True:
            existing = os.path.join(copydir, os.path.basename(x1d))
            if os.path.exists(existing):
                os.remove(existing)
        # else: normal copy error will be triggered below
        shutil.copy(x1d, copydir)
    print(f"Copied custom x1ds to {copydir}")


def apply_cos_shifts(infiledir, outdir, shift_file=None, targ=None, copydir=None, 
                     overwrite=False):
    if isinstance(infiledir, str):
        infiledir = []
    for item in infiledir:
        if os.path.isdir(item):
            if shift_file is None:
                targ = targ.lower()
                shift_file = determine_dir_shifts(item, targ)
            apply_shifts_dir(indir, outdir, shift_file)
        elif os.path.isfile(item):
            if shift_file is None:
                targ = targ.lower()
                shift_file = determine_file_shifts(item, targ)
            apply_shifts_file(infile, outdir, shift_file)
    add_hlsp_lvl0(outdir)
    if copydir is not None:
        copy_output_x1ds(outdir, copydir, overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="infiledir",
                        help="Input file or directory of files to shift")
    parser.add_argument(dest="outdir", 
                        help="Directory for output shifted 1D spectra")
    parser.add_argument("-s", "--shift_file", default=None,
                        help="Name of COS shifts file")
    parser.add_argument("-t", "--targ", default=None,
                        help="Name of target to shift")
    parser.add_argument("--copydir", default=None,
                        help="Name of directory to copy shifted products to")
    parser.add_argument("-c", "--overwrite", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")
    args = parser.parse_args()

    apply_cos_shifts(args.infiledir, args.outdir, args.shift_file, args.targ,
                     args.copydir, args.overwrite)
