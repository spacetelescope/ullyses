"""
To compare two HLSPs from different DRs:
python check_hlsps.py --oldhlsp <filename> --newhlsp <filename>
This prints the fitsdiff report and creates comparison plots

To compare a new HLSP to its contituent files:
python check_hlsps.py --newhlsp <filename>
This creates comparison plots only.
"""


import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import glob
import argparse

DRDIR = "/astro/ullyses/all_vetted_data_dr3"

def compare_hlsps(oldfile, newfile):
    # First do a fitsdiff
    fd = fits.FITSDiff(oldfile, newfile)
    print(fd.report())

    # Now make comparison plots
    old = fits.getdata(oldfile)
    new = fits.getdata(newfile)

    targ = fits.getval(newfile, "targname")

    # Plot wavelength vs flux, error, and SNR for both HLSPs
    for arr in ["flux", "error", "snr"]:
        plt.plot(old["wavelength"][0], old[arr][0], "#1b9e77", alpha=0.7, label="Old")
        plt.plot(new["wavelength"][0], new[arr][0], "#d95f02", alpha=0.7, label="New")
        plt.xlabel("Wavelength")
        plt.ylabel(arr)
        plt.title(targ)
        plt.legend()
        print("Close plot, then press enter to continue to next plot")
        plt.show()
        inp = input("")
        plt.clf()

def compare_x1ds(newhlsp):
    # Read in HLSP data and provenance tables
    hlsp = fits.getdata(newhlsp)
    prov = fits.getdata(newhlsp, 2)
    # These are the 1D spectra that contributed to the HLSP
    targ = fits.getval(newhlsp, "targname")
    x1d_filenames = prov['filename']
    x1ds = [os.path.join(DRDIR, targ.lower(), x1d) for x1d in x1d_filenames]

    # Plot wavelength vs flux and error for the HLSP and all constituent files
    for arr in ["flux", "error"]:
        plt.plot(hlsp["wavelength"][0], hlsp[arr][0], c="k", label="HLSP", zorder=100, alpha=0.6)
        for i,x1d in enumerate(x1ds):
            # Special handling for FUSE data
            if "vo.fits" in x1d:
                x1d = os.path.join(DRDIR, targ.lower(), f"dqscreened_{x1d_filenames[i].lower()}")
                sdq = 2
                wlarr = "wave"
                if arr == "error":
                    arr = "sigma"
            else:
                if arr == "sigma":
                    arr = "error"
                sdq = fits.getval(x1d, "sdqflags", 1)
                wlarr = "wavelength"
            
            # Plot x1d data
            x1ddata = fits.getdata(x1d)
            good_dq = np.where((x1ddata["dq"] & sdq) == 0)
            plt.plot(x1ddata[wlarr][good_dq], x1ddata[arr][good_dq], alpha=0.7, label=x1d_filenames[i])
        
        plt.xlabel("Wavelength")
        plt.ylabel(arr)
        plt.title(targ)
        plt.legend()
        print("Close plot, then press enter to continue to next plot")
        plt.show()
        inp = input("")
        plt.clf()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--oldhlsp", default=None,
                        help="Name of previous DR's HLSP")
    parser.add_argument("--newhlsp", default=None,
                        help="Name of upcoming DR's HLSP")
    args = parser.parse_args()

    oldhlsp = args.oldhlsp
    newhlsp = args.newhlsp
    assert newhlsp != None, "Must supply at least a new HLSP"
    if oldhlsp != None and newhlsp != None:
        compare_hlsps(oldhlsp, newhlsp)
    else:
        compare_x1ds(newhlsp)

