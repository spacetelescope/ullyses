import shutil
import glob
import os
import argparse
from astropy.io import fits

from ullyses_utils.parse_csv import parse_aliases
from ullyses_utils.fuse_target_info import FILESTOEDIT, FUSE_TARGS, TARGS_TO_FLAG
from fuse_add_dq import add_dq_col


"""
This script is used to flag bad wavelength regions
in FUSE data. The bad regions for each target are
listed in fuse_target_info.py in ullyses-utils.

Inputs:
    indir (str): Path to input directory
    outdir (str): Path to output directory
    overwrite (bool): If true, files in outdir will be overwritted. Default=False
"""


def flag_data(indir, outdir, overwrite=False):
    """
    This code steps through the targets and flags
    the DQs appropriately

    :param indir: Path to input directory
    :param outdir: Path to output directory
    :param overwrite: If true, overwrite output. Default=False
    :return:
    """

    targstoedit = list(FILESTOEDIT.keys())
    all_targs = []
    for k, v in FUSE_TARGS.items():
        all_targs += v
    for targ in all_targs:
        vofiles0 = glob.glob(os.path.join(indir, targ, "*_vo.fits"))
        vofiles = [x for x in vofiles0 if "dqscreened" not in x]
        if len(vofiles) > 1:
            print(f"more than one VO file for {targ}, exiting")
            break
        vofile = vofiles[0]
        vofilename = os.path.basename(vofile)
        outfilename = "dqscreened_" + vofilename
        outfile = os.path.join(outdir, targ, outfilename)
        # We don't edit FUSE data from DR to DR, so if a product already exists,
        # skip that target. This might change in the future.
        if os.path.exists(outfile):
            continue
        if targ in targstoedit:
            pars = FILESTOEDIT[targ]
            add_dq_col(vofile, outfile, pars["minwl"], pars["maxwl"], pars["dq"], overwrite=True)
        else:
            add_dq_col(vofile, outfile, [], [], [], overwrite=True)
    print("Flagged VO files")

    fuse_targname = fits.getval(vofile, "targname")
    aliases = parse_aliases()
    alias_mask = aliases.apply(lambda row: row.astype(str).str.fullmatch(re.escape(targ.upper())).any(), axis=1)
    if set(alias_mask) != {False}:
        ull_targname = aliases[alias_mask]["ULL_MAST_name"].values[0]
    else:
        raise KeyError(f"FUSE target {fuse_targname} not found in ULLYSES alias list")

    vofilename = os.path.basename(vofile)
    outfilename = "dqscreened_" + vofilename
    outfile = os.path.join(outdir, outfilename)
    if os.path.exists(outfile) and overwrite is True:
        os.remove(outfile)
    elif os.path.exists(outfile) and overwrite is False:
        print(f"Skipping {vofile}, output already exists and overwrite is False")
        return outfile

    targstoedit = list(TARGS_TO_FLAG.keys())
    if ull_targname in targstoedit:
        pars = TARGS_TO_FLAG[targ]
        add_dq_col(vofile, outfile, pars["minwl"], pars["maxwl"], pars["dq"], overwrite=overwrite)
    else:
        print(f"No bad regions in {vofile}, but still adding DQ array")
        add_dq_col(vofile, outfile, [], [], [], overwrite=True)


def copy_data(outdir, copydir):
    """
    Copies the output files into another directory.
    :param outdir: Path to the output directory
    :param copydir: Path to the directory to copy to
    :return: None
    """
    all_targs = []
    for k, v in FUSE_TARGS.items():
        all_targs += v
    for targ in all_targs:
        screened = glob.glob(os.path.join(outdir, targ, "dqscreened*_vo.fits"))
        for item in screened:
            destdir = os.path.join(copydir, targ.lower())
            if not os.path.exists(destdir):
                os.makedirs(destdir)
            shutil.copy(item, destdir)
    print(f"Copied flagged files to {destdir}")


def main(indir, outdir, overwrite):
    """
    Perform flagging.
    :param indir: Path to input directory
    :param outdir: Path to output directory
    :param overwrite: If true, overwrite output. Default=False
    :return: None
    """

    flag_data(indir, outdir, overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default=".",
                        help="Path to input directory with input NVO files to flag")
    parser.add_argument("-o", "--outdir",
                        help="Output path to place flagged NVO files")
    parser.add_argument("-c", "--overwrite", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")

    args = parser.parse_args()

    main(args.indir, args.outdir, args.overwrite)
