import argparse

from ullyses.create_ullyses_hlsp import make_lcogt_tss


def lcogt_hlsps_wrapper(indir, outdir, targ, hlspname=None, photfile=None):
    make_lcogt_tss(indir=indir, outdir=outdir, targ=targ, photfile=photfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",
                        help="Name of directory that houses FITS files used to perform photometry")
    parser.add_argument("-o", "--outdir",
                        help="Directory to write HLSPs to")
    parser.add_argument("-t", "--targ",
                        help="ULLYSES target name")
    parser.add_argument("--hlspname", default=None,
                        help="Name of output HLSP file. By default, follows ULLYSES standard")
    parser.add_argument("--photfile", default=None,
                        help="Name of file containing photometric measurements")
    args = parser.parse_args()

    make_lcogt_tss(args.indir, args.outdir, args.targ, args.hlspname, args.photfile)
