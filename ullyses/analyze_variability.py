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

COLORS = ["#9970ab", "#5aae61", "#d95f02", "#e7298a", "#1bddf2", "#1f78b4", "#fb9a99",
          "#fdbf6f", "#ffe44b", "#b15928", "#cab2d6", "#b2df8a", "#749688", "#7a7a7a",
          "#911cff", "#a6cee3"]

#-------------------------------------------------------------------------------

def add_tts_regions(fig, detector, wmin, wmax):
    ''' Adding highlights to regions of interest for the T-Tauri stars.
    '''

    if 'FUV' in detector:
        #from Section 1 of Ardila et al. 2013, ApJS, 207, 1
        regions = [(1238, 1243), # N V doublet: 1239 + 1243
                   (1393, 1404), # Si IV doublet: 1394 + 1403
                   (1547, 1552), # C IV doublet: 1548 + 1551
                   (1638, 1642) # He II 1640
                   ]
    elif 'NUV' in detector:
        #from Table 5 of Robinson & Espaillat 2019, ApJ, 874, 129
        regions = [(2323, 2327), # C II 2325
                   (2668, 2672), # Al III 2670
                   (2794, 2798), # Mg II 2796
                   ]
    elif detector == 'CCD':
        # from Table B.1 of Alcala et al. 2017, A&A, 600, A20
        regions = [(4860, 4864), # Hβ 4862
                   (6561, 6565), # Hα 6563
                   (5874, 5878), # He I 5876
                   (7771, 7775), #O I 7773
                   ]
    else:
        print(f'Unrecognized Detector: {detector}')
        regions = []

    # plot the regions
    for x0, x1 in regions:
        if x0 > wmax or x1 < wmin:
            # don't plot regions outside of the wavelength range
            continue
        fig.add_vrect(x0=x0, x1=x1, line_width=0, fillcolor="mistyrose", opacity=0.2, layer='below')

    return fig

#-------------------------------------------------------------------------------

def compare_spectra(files, use_grating=None, savefig=True, savedir="", tts_regions=True):

    # sort the files by grating
    files_d = defaultdict(list)
    for item in files:
        files_d[fits.getval(item, "opt_elem")].append(item)

    file_targets = np.unique([match_aliases(fits.getval(f, 'TARGNAME')) for f in files])
    if len(file_targets) > 1:
        raise ValueError(f'Too many targets in filelist. Expecting only 1: {file_targets}')
    elif len(file_targets) == 0:
        print('No Files!')
        return
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

        median_flux = [] # used for calculating the median value of all files
        median_wave = [] # used for plotting the median wave value of all files

        # plot the files
        for i, item in enumerate(grating_files):
            n_ext = fits.getval(item, "nextend") + 1 # how many extensions are in the file
            date_obs = fits.getval(item, "date-obs", 1) # date of observation

            # open the file to plot
            with fits.open(item) as hdulist:
                detector = hdulist[0].header['DETECTOR']
                for ext in range(1, n_ext):
                    # make sure to plot every SCI extension!
                    if hdulist[ext].name != "SCI":
                        continue
                    data = hdulist[ext].data

                    wl = data["wavelength"].flatten()
                    flux = data["flux"].flatten()

                    # add the flux to the median array to be calculated
                    median_flux.append(flux)
                    median_wave.append(wl)

                    lbl = f"{os.path.basename(item)} {date_obs} ext{ext}"

                    fig.add_trace(go.Scatter(x=wl,
                                             y=flux,
                                             mode='lines',
                                             line=dict(color=COLORS[i]),
                                             opacity=0.8,
                                             name=lbl,
                                             ))
        if tts_regions:
            # plot boxes around regions of interest for TTS
            fig = add_tts_regions(fig, detector, wl.min(), wl.max())

        # calculate the median of all files to overplot
        fig.add_trace(go.Scatter(x=np.median(median_wave, axis=0),
                                 y=np.median(median_flux, axis=0),
                                 mode='lines',
                                 line=dict(color='gainsboro',
                                           width=6),
                                 opacity=0.8,
                                 name='Median Value'
                                 ) )

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
            os.chmod(figname, 0o774)
            print(f"Saved {figname}")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", type=str, default='/astro/ullyses/stis_data_checks/',
                        help="Top level directory to search for files")
    parser.add_argument("-o", "--outdir", type=str, default='',
                        help="Directory to save the outputs. Default behavior is to save to the target input directory")
    parser.add_argument("-t", "--target", type=str, default="",
                        help="Use if you are specifying a single target directory")
    parser.add_argument("-g", "--grating", type=str, default = None,
                        help="Specify a grating")
    parser.add_argument("-s", "--no_save", action='store_false', default=True,
                        help="use if you do not want to save a plot")
    args = parser.parse_args()

    if args.target != "":
        targ_dirs = [os.path.join(args.indir, args.target.lower())]
    else:
        targ_dirs = glob.glob(os.path.join(args.indir, '*'))

    for targ_dir in targ_dirs:
        print(targ_dir)
        if not os.path.isdir(targ_dir):
            # skipping over any files that live in the top level directory
            continue

        files = np.sort(glob.glob(os.path.join(targ_dir, '*x1d*')))

        if args.outdir == '':
            # if no outdir is specified, save them in the target directory
            compare_spectra(files, use_grating=args.grating, savefig=args.no_save, savedir=targ_dir)
        else:
            # otherwise, save the files to the specified directory
            if not os.path.exists(args.outdir):
                # make the directory if it doesn't exist
                os.mkdir(args.outdir)
                print(f'Created directory: {args.outdir}')
            compare_spectra(files, use_grating=args.grating, savefig=args.no_save, savedir=args.outdir)
