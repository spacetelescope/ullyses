"""
Overplot multiple spectra for a single target, to see if there is flux
variability with time. If there is, the ULLYSES team will make a time-series
product for the target.
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
plt.style.use('tableau-colorblind10')
import os
import glob
import argparse
from collections import defaultdict

COLORS = ["#9970ab", "#5aae61", "#d95f02", "#e7298a", "#1bddf2", "#1f78b4", "#fb9a99", "#fdbf6f", "#ffe44b", 
          "#b15928", "#cab2d6", "#b2df8a", "#000000", "#7a7a7a", "#911cff", "#a6cee3"]


def compare_spectra(files, targname="", use_grating=None, savefig=True):
    fig,ax = plt.subplots(figsize=(12,6))
    files_d = defaultdict(list)
    for item in files:
        files_d[fits.getval(item, "opt_elem")].append(item)
    for grating,files in files_d.items():
        if use_grating is not None and grating != use_grating:
            continue
        for i,item in enumerate(files):
            n_ext = fits.getval(item, "nextend") + 1
            date_obs = fits.getval(item, "date-obs", 1)
            with fits.open(item) as hdulist:
                for ext in range(1, n_ext):
                    if hdulist[ext].name != "SCI":
                        continue
                    data = hdulist[ext].data
                    wl = data["wavelength"].flatten()
                    flux = data["flux"].flatten()
                    lbl = f"{os.path.basename(item)} {date_obs}"
                    ax.plot(wl, flux, color=COLORS[i], alpha=.8, label=lbl)
        ax.set_xlabel("Wavelength [A]")
        ax.set_ylabel("Flux [erg/s/cm**2/Angstrom]")
        ax.legend()
        title = f"{targname} {grating}"
        ax.set_title(title)
        if savefig is True:
            figname = f"{targname}_{grating}.png"
            plt.savefig(figname)
            print(f"Saved {figname}")

