"""
Overplot multiple spectra for a single target, to see if there is flux
variability with time. If there is, the ULLYSES team will make a time-series
product for the target.
"""

import argparse
import glob
import os
import numpy as np
import plotly.graph_objects as go

from astropy.io import fits
from collections import defaultdict

from ullyses_utils.match_aliases import match_aliases

COLORS = ["#9970ab", "#5aae61", "#d95f02", "#e7298a", "#1bddf2", "#1f78b4", "#fb9a99", "#fdbf6f", "#ffe44b",
          "#b15928", "#cab2d6", "#b2df8a", "#000000", "#7a7a7a", "#911cff", "#a6cee3"]

def compare_spectra(files, use_grating=None, savefig=True, savedir=""):

    # sort the files by grating
    files_d = defaultdict(list)
    for item in files:
        files_d[fits.getval(item, "opt_elem")].append(item)

    file_targets = np.unique([match_aliases(fits.getval(f, 'TARGNAME')) for f in files])
    if len(file_targets) > 1:
        raise ValueError(f'Too many targets in filelist. Expecting only 1: {file_tagets}')
    elif len(file_targets) == 0:
        raise ValueError(f'No files!')
    else:
        targname = file_targets[0]

    for grating, grating_files in files_d.items():

        # initialize figure
        fig = go.Figure()

        if use_grating is not None and grating != use_grating:
            # skip the file for now if the grating is not used
            print(f"{grating} does not match expected {use_grating}. Skipping.")
            continue

        if len(grating_files) <= 1:
            # we don't need to check the variability if there is only 1 file
            print(f'Only 1 file: {grating_files}. Skipping {grating}.')
            continue

        print(f'Plotting {len(grating_files)} files for {targname} {grating}')

        # make sure there are enough colors to plot with
        #   this should not be needed very often
        global COLORS
        while len(COLORS) < len(grating_files):
            COLORS = COLORS * 2

        # plot the files
        for i, item in enumerate(grating_files):
            n_ext = fits.getval(item, "nextend") + 1 # how many extensions are in the file
            date_obs = fits.getval(item, "date-obs", 1) # date of observation

            # open the file to plot
            with fits.open(item) as hdulist:
                for ext in range(1, n_ext):
                    # make sure to plot every SCI extension!
                    if hdulist[ext].name != "SCI":
                        continue
                    data = hdulist[ext].data

                    wl = data["wavelength"].flatten()
                    flux = data["flux"].flatten()

                    lbl = f"{os.path.basename(item)} {date_obs} ext{ext}"

                    fig.add_trace(go.Scatter(x=wl,
                                             y=flux,
                                             mode='lines',
                                             line=dict(color=COLORS[i]),
                                             opacity=0.8,
                                             name=lbl,
                                             ))

        # set up plot labels
        fig.update_layout(yaxis={'showexponent': 'all',
                                 'exponentformat': 'E',
                                 'rangemode': 'tozero'},
                          yaxis_title={'text' : "Flux [erg/s/cm<sup>2</sup>/&#8491;]",
                                       'font' : {'size' : 20}},
                          xaxis_title={'text' : "Wavelength [&#8491;]",
                                       'font' : {'size' : 20}},
                          title={'text' : f"{targname} {grating}",
                                 'font' : {'size' : 30}},
                          template='plotly_dark'
                          )

        if savefig is True:
            figname = os.path.join(savedir, f"{targname}_{grating}_timeseries_check.html")
            fig.write_html(figname)
            os.chmod(filename, 774)
            print(f"Saved {figname}")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dir", type=str, default='/astro/ullyses/stis_data_checks/',
                        help="Top level directory to search for files")
    parser.add_argument("-t", "--target", type=str, default="",
                        help="Use if you are specifying a single target directory")
    parser.add_argument("-g", "--grating", type=str, default = None,
                        help="Specify a grating")
    parser.add_argument("-s", "--no_save", action='store_false', default=True,
                        help="use if you do not want to save a plot")
    args = parser.parse_args()

    if args.target != "":
        targ_dirs = [os.path.join(args.dir, args.target)]
    else:
        targ_dirs = glob.glob(os.path.join(args.dir, '*'))

    for targ_dir in targ_dirs:
        print(targ_dir)
        if not os.path.isdir(targ_dir):
            # skipping over any files that live in the top level directory
            continue

        files = glob.glob(os.path.join(targ_dir, '*x1d*'))

        compare_spectra(files, use_grating=args.grating, savefig=args.no_save, savedir=targ_dir)
