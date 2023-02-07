"""
COS NUV data are vignetted in the beginning of all channels. For monitoring stars,
we actually correct the vignetting because it encroaches on the MgII lines in 
G230L/2950 using overlapping G230L/2635 data.
But for archival data, we do not necessarily have overlapping data to correct
vignetting, and also we do not necessarily care about vignetting effects if there
are no lines in the vignetted region.
This code can be used to flag the DQ array in the vignetted region in COS/NUV data. 
"""

from astropy.io import fits
import argparse


def flag_vignetting(filename, flag_npixels=200):
    """Flag vignetted region of each COS/NUV stripe as DQ=8 (poorly calibrated region).
    DQ=8 is in SDQFLAGS.

    Args:
        filename (str): Name of file to be updated.
        flag_npixels (int): Number of pixels on the left edge of each stripe to flag.
             Default is 200. 
    """

    hdr0 = fits.getheader(filename, 0)
    assert (hdr0["instrume"] == "COS" and hdr0["detector"] == "NUV"), f"Input filename is not COS/NUV: {filename}"
    with fits.open(filename, mode="update") as hdulist:
        hdulist[1].data["DQ"][0][:flag_npixels] |= 8
        hdulist[1].data["DQ"][1][:flag_npixels] |= 8
        hdulist[1].data["DQ"][2][:flag_npixels] |= 8
    print(f"DQ=8 added to first {flag_npixels} pixels for each strip in file {filename}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="filename",
                        help="Input file to flag")
    parser.add_argument("-n", "--flag_npixels", default=200,
                        help="Number of pixels to flag on each stripe")
    args = parser.parse_args()

    flag_vignetting(args.filename, args.flag_npixels)

