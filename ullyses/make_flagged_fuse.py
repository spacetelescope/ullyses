"""
This script is used to flag bad wavelength regions
in FUSE data. The bad regions for each target are
listed in TARGS_TO_FLAG.


Inputs:
    indir (str): Path to input directory
    outdir (str): Path to output directory
    overwrite (bool): If true, files in outdir will be overwritted. Default=False
"""

import re
import shutil
import glob
import os
import argparse
from astropy.io import fits

from ullyses.fuse_add_dq import add_dq_col
import ullyses_utils
from ullyses_utils.parse_csv import parse_aliases
from ullyses_utils.readwrite_yaml import read_config
utils_dir = ullyses_utils.__path__[0]
fuse_dir = os.path.join(utils_dir, "data", "fuse")
TARGS_TO_FLAG = read_config(os.path.join(fuse_dir, "fuse_dq_flagging.yaml"))


def flag_file(vofile, outdir, ull_targname=None, overwrite=False):
    """
    This code steps through the targets and flags
    the DQs appropriately

    :param vofile: Path to VO file
    :param outdir: Path to output directory
    :param ull_targname: If known, name of official ULLYSES MAST target name.
    :param overwrite: If true, overwrite output. Default=False
    :return: None
    """

    if "dqscreened_" in vofile:
        print(f"Input file {vofile} is already flagged, skipping")
        return vofile

    if ull_targname is None:
        fuse_targname = fits.getval(vofile, "targname")
        aliases = parse_aliases()
        alias_mask = aliases.apply(lambda row: row.astype(str).str.fullmatch(re.escape(fuse_targname.upper())).any(), axis=1)
        if set(alias_mask) != {False}:
            ull_targname = aliases[alias_mask]["target_name_hlsp"].values[0]
        else:
            raise KeyError(f"FUSE target {fuse_targname} not found in ULLYSES alias list")

    if not os.path.exists(outdir):
        os.makedirs(outdir)
    vofilename = os.path.basename(vofile)
    outfilename = "dqscreened_" + vofilename
    if outfilename.endswith(".fit"):
        outfilename = outfilename.replace(".fit", ".fits")
    outfile = os.path.join(outdir, outfilename)
    if os.path.exists(outfile) and overwrite is True:
        os.remove(outfile)
    elif os.path.exists(outfile) and overwrite is False:
        print(f"Skipping {vofile}, output already exists and overwrite is False")
        return outfile

    targstoedit = list(TARGS_TO_FLAG.keys())
    if ull_targname in targstoedit:
        pars = TARGS_TO_FLAG[ull_targname]
        add_dq_col(vofile, outfile, pars["minwl"], pars["maxwl"], pars["dq"], overwrite=overwrite)
    else:
        print(f"No bad regions in {vofile}, but still adding DQ array")
        add_dq_col(vofile, outfile, [], [], [], overwrite=True)

    return outfile


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", default=".",
                        help="Path to input directory with input NVO files to flag")
    parser.add_argument("-o", "--outdir",
                        help="Output path to place flagged NVO files")
    parser.add_argument("-t", "--targ",
                        help="ULLYSES name of target")
    parser.add_argument("-c", "--overwrite", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")

    args = parser.parse_args()

    flag_file(args.infile, args.outdir, args.targ, args.overwrite)
