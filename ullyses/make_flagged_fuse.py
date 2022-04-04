import shutil
import glob
import os
import argparse
from ullyses_utils.fuse_target_info import FILESTOEDIT, FUSE_TARGS
from fuse_add_dq import add_dq_col


"""
This script is used to flag bad wavelength regions
in FUSE data. The bad regions for each target are
listed in fuse_target_info.py in ullyses-utils.

Inputs:
    indir: (str) Path to input directory
    outdir: (str) Path to output directory
    copydir: (str) Path to directory to copy output to, if specified
"""


def flag_data(indir, outdir):
    """
    This code steps through the targets and flags
    the DQs appropriately


    :param indir:
    :param outdir:
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


def main(indir, outdir, copydir):
    """
    Perform flagging.
    :param indir: Path to input directory
    :param outdir: Path to output directory
    :param copydir: Path to directory to copy output to, if specified
    :return: None
    """

    flag_data(indir, outdir)
    if copydir is not None:
        copy_data(outdir, copydir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default=".",
                        help="Path to input directory with ASN and rawtag files")
    parser.add_argument("-o", "--outdir",
                        help="Path to output directory")
    parser.add_argument("-c", "--copydir", default=None,
                        help="Path to directory to copy x1ds to, if so desired.")

    args = parser.parse_args()

    main(args.indir, args.outdir, args.copydir)
