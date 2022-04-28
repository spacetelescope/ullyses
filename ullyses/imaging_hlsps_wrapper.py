import argparse

from ullyses.create_ullyses_hlsp import make_imaging_hlsps

def imaging_hlsps_wrapper(drcfile, outdir, targ, hdr_targ=None, hlspname=None):
    make_imaging_hlsps(drcfile=item, outdir=outdir, targ=targ, hdr_targ="NGC3109")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="drcfile",
                        help="Name of DRC file to turn into an HLSP")
    parser.add_argument("-o", "--outdir",
                        help="Directory to write HLSP to")
    parser.add_argument("-t", "--targ",
                        help="ULLYSES target name")
    parser.add_argument("--hdr_targ", default=None,
                        help="If specified, alternative target name to use in HLSP file name")
    parser.add_argument("--hlspname", default=None,
                        help="Name of output HLSP file. By default, follows ULLYSES standard")
    args = parser.parse_args()

    imaging_hlsps_wrapper(args.drcfile, args.outdir, args.targ, args.hdr_targ, args.hlspname)
