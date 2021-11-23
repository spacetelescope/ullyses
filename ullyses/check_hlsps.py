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

DRDIR = "/astro/ullyses/all_vetted_data_dr4"

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
    x1d_rootnames = [x.split("_")[0] for x in x1d_filenames]
    x1ds = []
    for item in x1d_rootnames:
        x1d = glob.glob(os.path.join(DRDIR, targ.lower(), "*"+item+"*.fits"))
        assert len(x1d) == 1, f"Expected one matching file for {item}, got {len(x1d)}"
        x1ds.append(x1d[0])

    # Plot wavelength vs flux and error for the HLSP and all constituent files
    for arr in ["flux", "error"]:
        plt.plot(hlsp["wavelength"][0], hlsp[arr][0], c="k", label="HLSP", zorder=100, alpha=0.6)
        for i,x1d in enumerate(x1ds):
            # Special handling for FUSE data
            if "vo.fits" in x1d:
                x1d = os.path.join(DRDIR, targ.lower(), f"dqscreened_{x1d_filenames[i].lower()}")
                sdq = 3
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


def compare_tss(newtss, oldtss=None):
    newdata = fits.getdata(newtss)
    newwl = newdata["wavelength"][0]
    if oldtss is not None:
        olddata = fits.getdata(oldtss)
        oldwl = olddata["wavelength"][0]
        assert len(newdata["flux"][0]) == len(olddata["flux"][0]), \
            f"{newtss}: New and old TSS dimensions do not match!"
    l = len(newdata["flux"][0])
    inds = np.arange(0, l, l/5, dtype=int)
    for i in inds:
        newflux = newdata["flux"][0][i]
        newerr = newdata["error"][0][i]
        fig, axes0 = plt.subplots(2, 1, figsize=(15,10), sharex=True)
        axes = axes0.flatten()
        axes[0].plot(newwl, newflux, "r", alpha=0.8, label="New")
        axes[0].set_title("Flux")
        axes[0].set_ylabel("Flux")
        if oldtss is not None:
            oldflux = olddata["flux"][0][i]
            olderr = olddata["error"][0][i]
            oldsnr = oldflux / olderr
            axes[0].plot(oldwl, oldflux, "k", label="Old")
            axes[1].plot(oldwl, oldsnr, "k", label="Old")
            axes[0].legend(loc="upper right")
            axes[1].legend(loc="upper right")
        newsnr = newflux / newerr
        axes[1].plot(newwl, newsnr, "r", alpha=0.8, label="New")
        axes[1].set_title("SNR")
        axes[1].set_xlabel("Wavelength")
        axes[1].set_ylabel("SNR")
         
        mjdstart = newdata["mjdstart"][0][i]
        mjdend = newdata["mjdend"][0][i]
        fig.suptitle(f"Index={i}, {mjdstart:.4f} - {mjdend:.4f}")
        print("Close plot, then press enter to continue to next plot")
        plt.show()
        inp = input("")


def check_tss(newhlsp):
    newdata = fits.getdata(newhlsp)
    opt_elem = fits.getval(newhlsp, "dispersr")
    newwl = newdata["wavelength"][0]
    l = len(newdata["flux"][0])
    fig, axes0 = plt.subplots(2, 1, figsize=(15,10), sharex=True)
    for i in range(l):
        newflux = newdata["flux"][0][i]
        newerr = newdata["error"][0][i]
        newsnr = newflux / newerr
        axes = axes0.flatten()
        axes[0].plot(newwl, newflux, alpha=0.8, label=f"Index={i}")
        axes[0].set_title("Flux")
        axes[0].set_xlabel("Wavelenth")
        axes[0].set_ylabel("Flux")
        if opt_elem == "G230L":
            axes[0].set_xlim(2790, 2810)
            axes[0].set_ylim(top=1.4e-12)
        else:
            axes[0].set_xlim(1545, 1557)

        axes[1].plot(newwl, newsnr, alpha=0.8, label=f"Index={i}")
        axes[1].set_title("SNR")
        axes[1].set_xlabel("Wavelenth")
        axes[1].set_ylabel("SNR")
        if i % 10 == 0 and i != 0:
            axes[0].legend(loc="upper right")
            axes[1].legend(loc="upper right")
            print("Close plot, then press enter to continue to next plot")
            plt.show()
            inp = input("")
            fig, axes0 = plt.subplots(2, 1, figsize=(15,10), sharex=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--oldhlsp", default=None,
                        help="Name of previous DR's HLSP")
    parser.add_argument("--newhlsp", default=None,
                        help="Name of upcoming DR's HLSP")
    parser.add_argument("--tss", default=False, action="store_true",
                        help="If toggled, the input HLSPs are time-series spectra")
    args = parser.parse_args()

    oldhlsp = args.oldhlsp
    newhlsp = args.newhlsp
    tss = args.tss
    assert newhlsp != None, "Must supply at least a new HLSP"
    if oldhlsp != None and newhlsp != None and tss is False:
        compare_hlsps(oldhlsp, newhlsp)
    if oldhlsp != None and newhlsp != None and tss is True:
        compare_tss(newhlsp, oldhlsp) 
    elif tss is True:
        check_tss(newhlsp)
    else:
        compare_x1ds(newhlsp)

