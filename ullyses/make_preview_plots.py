import argparse
import glob
import os
import numpy as np
import plotly.graph_objects as go
import astropy.units as u
from plotly.subplots import make_subplots
from itertools import repeat

from astropy.io import fits
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region

from ullyses_utils.match_aliases import match_aliases
from ullyses_utils.parse_csv import parse_database_csv

#-------------------------------------------------------------------------------

def create_visibility(trace_lengths, visible_list):
    '''Create visibility lists for plotly buttons. trace_lengths and visible_list
    must be in the correct order.

    Inputs:
    -------
    - trace_lengths : list
        List of the number of traces in each "button set".
    - visible_list : list
        Visibility setting for each button set (either True or False).

    Outputs:
    --------
    - visibility : list
        an updated list of the visibility for each added trace
    '''

    visibility = []  # Total visibility. Length should match the total number of traces in the figure.
    for visible, trace_length in zip(visible_list, trace_lengths):
        visibility += list(repeat(visible, trace_length))  # Set each trace per button.

    return visibility

#-------------------------------------------------------------------------------

def make_buttons(fig, all_vis, n_all_traces, button_names):
    '''Create buttons for plotly figures, being careful to count all of the
    number of traces being added.

    Inputs:
    -------
    - fig : plotly figure
        the figure with all traces added
    - all_vis : list
        A list of the trace visibilities for which traces should shown up
        when a button is selected
    - n_all_traces : int
        A count of the number of added traces
    - button_names : list
        What should be read on the buttons

    Outputs:
    --------
    - fig : plotly figure
        the figure returned with updated buttons
    '''

    # Create a dictionary of the visibility parameters based on the number of traces
    # all of the traces that are drawn per subplot plus what we want to show by defualt
    # vis will look like [True, False, False], [False, True, False], etc
    all_vis = [create_visibility(n_all_traces, vis) for vis in all_vis]

    buttons = [dict(label=str(label),
                    method='relayout',
                    args=[{'visible': visibility,
                           'yaxis.type' : scale_type,
                           'yaxis2.type' : scale_type,
                           },])
               for label, visibility, scale_type in zip(button_names, all_vis, ['linear', 'log'])]

    # Update remaining layout properties
    fig.update_layout(updatemenus=[dict(type="buttons",
                                        direction="right",
                                        x=0.98,
                                        y=1.05,
                                        showactive=True,
                                        buttons=buttons,
                                        font=dict(color='black',
                                                  size=16),
                                    )])

    return fig

#-------------------------------------------------------------------------------

def add_tts_regions(fig, detector):
    ''' Adding highlights to regions of interest for the T-Tauri stars.
    '''

    x_mins = []
    x_maxs = []
    for trace_data in fig.data:
        x_mins.append(min(trace_data.x))
        x_maxs.append(max(trace_data.x))
    wmin = min(x_mins)
    wmax = max(x_maxs)

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
                   (7771, 7775), # O I 7773
                   ]
    else:
        print(f'Unrecognized Detector: {detector}')
        regions = []

    # plot the regions
    for x0, x1 in regions:
        if x0 > wmax or x1 < wmin:
            # don't plot regions outside of the wavelength range
            continue
        fig.add_vrect(x0=x0, x1=x1, line_width=0, fillcolor="black", opacity=0.2, layer='below')

    return fig

#-------------------------------------------------------------------------------

def get_max_flux(wave, flux, dlam=10., echelle=False):
    '''
    Determine the maximum flux in spectrum after excluding the contribution
      from strong airglow lines. Currently the absolute maximum is chosen,
      so the value is not robust against noise.
    '''

    airglow_lines = {949.7430 : 'H I Lyman Delta',
                     972.5367 : 'H I Lyman Gamma',
                     988.773 : 'O I UV 5',
                     1025.7219 : 'H I Lyman Beta',
                     1168.000 : 'He I UV 1 in 2nd order',
                     1215.6683 : 'H I Lyman Alpha',
                     1302.1685 : 'O I UV2',
                     1355.5977 : 'O I UV1',
                    }

    for i, line in enumerate(airglow_lines.keys()):
        # find where the airglow is in the spectrum (if any)
        #  Full width is dlam, so range is (lam0-dlam/2, lam0+dlam/2)
        wh_airglow = np.where( (wave >= line-(dlam/2)) & (wave <= line+(dlam/2)) )[0]
        if len(wh_airglow):
            # set the flux to zero in that region
            flux[wh_airglow] = 0.

    # exclude the edges as well;
    if echelle:
        # the echelle edges extend pretty far
        flux = flux[3000:-1001]
    elif len(flux) > 3000:
        flux = flux[1000:-1001]
    else:
        flux = flux[100:-100]

    return flux.max()

#-------------------------------------------------------------------------------

def add_spectral_trace(fig, wave, flux, lbl, row, legendgroup, color, showlegend=True, dash=None, visible=True):

    fig.add_trace(go.Scattergl(x=wave,
                               y=flux,
                               mode='lines',
                               line=dict(dash=dash,
                                         color=color,
                                         width=3),
                               opacity=0.6,
                               name=lbl,
                               showlegend=showlegend,
                               legendgroup=legendgroup,
                               visible=visible,
                               hovertemplate = '%{text}<br>(%{x}, %{y})<extra></extra>',
                               text = [f"<b>{'/'.join(lbl.split(' '))}</b>"]*len(wave),
                               ), row=row, col=1)

    return fig

#-------------------------------------------------------------------------------

def open_files_and_plot(fig, cspec_file, row, flux_maxes, legendgroup, color, visible=True, dash=None):

    #print(cspec_file)
    trace_count = 0
    with fits.open(cspec_file) as hdu:
        detector = hdu[0].header['DETECTOR']
        echelle = any([g.startswith('E') for g in hdu[2].data['DISPERSER']])
        all_gratings = [g for g in hdu[2].data['DISPERSER']]
        lbl = f"{hdu[0].header['TELESCOP']} {hdu[0].header['INSTRUME']} {'&'.join(np.unique(all_gratings))}"

        for ext in range(len(hdu[1].data)):
            wave = hdu[1].data["WAVELENGTH"][ext]
            flux = hdu[1].data["FLUX"][ext]

            if ext > 1:
                showlegend = False
            else:
                showlegend = True

            flux_maxes.append(get_max_flux(wave, flux, echelle=echelle))

            fig = add_spectral_trace(fig, wave, flux, lbl, row, legendgroup, color,
                                     showlegend=showlegend, dash=dash, visible=visible)
            trace_count += 1

    return fig, flux_maxes, detector, trace_count

#-------------------------------------------------------------------------------

def plot_preview(fig, targ_dir, visible=True):

    ntraces = 0

    prev = glob.glob(os.path.join(targ_dir, '*_preview-spec.fits'))
    if len(prev) == 0:
        return fig, ntraces

    with fits.open(prev[0]) as hdu:
        wave = hdu[1].data["WAVELENGTH"].ravel()
        flux = hdu[1].data["FLUX"].ravel()
        #flux[np.where(flux == 0.0)[0]] = "nan"

    for row in [1, 2]:

        if row == 1:
            showlegend = True
        else:
            showlegend = False

        fig.add_trace(go.Scatter(x=wave,
                                 y=flux,
                                 mode='lines',
                                 line=dict(width=3,
                                           color='dimgrey'),
                                 opacity=0.7,
                                 name='All Abutted Spectra<br>(Level 4 HLSP)',
                                 legendgroup='preview',
                                 showlegend=showlegend,
                                 visible=visible,
                                 hovertemplate = '<b>All Abutted</b><br>(%{x}, %{y})<extra></extra>',
                                 legendrank=1001, # in legend top if legend rank < 1000 & bottom if rank > 1000
                                 ), row=row, col=1)

        ntraces += 1

    return fig, ntraces

#-------------------------------------------------------------------------------

def update_plot_layouts(fig, ymax_arr, target):

    fig.update_annotations(font_size=25)

    fig.add_annotation(text="<sup>*</sup>Gratings with similar resolutions ",
                       xref="paper", yref="paper",
                       x=0.97, y=0.45, showarrow=False, font_size=18)

    fig.update_layout(plot_bgcolor='gainsboro',
                      yaxis={'showexponent': 'all',
                             'exponentformat': 'E',
                             'autorange' : False,
                             'rangemode': 'tozero',
                             'range' : [-0.1*np.max(ymax_arr), 1.10*np.max(ymax_arr)],
                             'title': {'text' : 'Flux (erg/s/cm<sup>2</sup>/&#8491;)',
                                       'font' : {'size' : 25}, },
                             'tickfont' : {'size' : 20},
                             },
                      yaxis2={'showexponent': 'all',
                              'exponentformat': 'E',
                              'autorange' : False,
                              'rangemode': 'tozero',
                              'range' : [-0.1*np.max(ymax_arr), 1.10*np.max(ymax_arr)],
                              'title': {'text' : 'Flux (erg/s/cm<sup>2</sup>/&#8491;)',
                                        'font' : {'size' : 25}, },
                              'tickfont' : {'size' : 20},
                              },
                      xaxis2={'title': {'text' : 'Wavelength (&#8491;)',
                                        'font' : {'size' : 25}, },
                              'tickfont' : {'size' : 20},
                              },
                      title={'text' : f"ULLYSES Products Preview for {target}",
                             'font' : {'size' : 30}},
                      hovermode='x unified',
                      )

    return fig

#-------------------------------------------------------------------------------

def high_and_low_df():

    high_csv, highmass_df = parse_database_csv('lmc')
    low_csv, lowmass_df = parse_database_csv('tts')

    for df in [highmass_df[0], lowmass_df[0]]:
        temp_list = []
        for t in df['target_name_ullyses']:
            temp_list.append(match_aliases(t))
        df['target_name_hlsp'] = temp_list

    return highmass_df[0], lowmass_df[0]

#-------------------------------------------------------------------------------

def make_plots(outdir_name, dr):

    # There are different levels for different combinations of gratings:
    #   panel 1 is individual gratings
    #   panel 2 is combined gratings
    #   overplotted on both is the "preview_spec" which can be turned on&off

    highmass_df, lowmass_df = high_and_low_df()

    # v-cv-cha
    # hd-104237e
    # av-456
    no_previews = []
    for targ_dir in np.sort(glob.glob("/astro/ullyses/ULLYSES_HLSP/*")):

        targ_dir = os.path.join(targ_dir, "dr7")

        ## sort out the different types of cspec files
        # CSPEC files: *_<target>_<grating>_dr?_cspec.fits
        # ASPEC files: *_<target>_<grating>_dr?_cspec.fits
        grating_cspec = glob.glob(os.path.join(targ_dir, '*_cspec.fits'))
        combined_aspec = glob.glob(os.path.join(targ_dir, '*_aspec.fits'))

        # check directories where there are no cspec files (should be tss only)
        if not len(grating_cspec):
            no_previews.append(targ_dir)
            continue

        # initialize a file
        fig = make_subplots(rows=2,
                            cols=1,
                            shared_xaxes=True,
                            #shared_yaxes=True, # this only shares across columns...
                            subplot_titles=('Single Grating Coaddition',
                                            'Single Instrument, Multiple Grating<sup>*</sup> Abutment'),
                            row_titles=('<a href="https://ullyses.stsci.edu/ullyses-data-description.html#productDescrip">Level 2 HLSP</a>',
                                        '<a href="https://ullyses.stsci.edu/ullyses-data-description.html#productDescrip">Level 3 HLSP</a>'),
                            vertical_spacing=0.05)

        # get the ULLYSES HLSP name for the target
        target = match_aliases(fits.getval(grating_cspec[0], 'TARGNAME'))

        # create empty buttons for multiple scale options
        buttons = ['Linear Scale', 'Log Scale']
        n_all_traces = []
        all_vis = []

        # looping through the traces to make my buttons. Visible by default for linear
        #   first. Not visible for log scale
        for i, (buttonname, vis) in enumerate(zip(buttons, [True, False])):
            # initializing a False array except for the button that will be visible
            #    will be adjusted to account for the number of traces per plot later
            current_plot_visibility = [False] * len(buttons)
            current_plot_visibility[i] = True
            all_vis.append(current_plot_visibility)

            temp_trace_count = 0

            ## row 1; grating only files (cspec)
            grating_ymaxes = []
            detectors = []
            for i_cspec, cspec_file in enumerate(np.sort(grating_cspec)):
                fig, grating_ymaxes, d, tc = open_files_and_plot(fig, cspec_file, 1, grating_ymaxes, f'{i_cspec}grating', CSPEC_COLORS[i_cspec], visible=vis)
                detectors.append(d)
                temp_trace_count += tc

            ## row 2; abutted files
            combined_maxes = []
            for i_aspec, aspec_file in enumerate(np.sort(combined_aspec)):
                fig, combined_maxes, d, tc = open_files_and_plot(fig, aspec_file, 2, combined_maxes, f'{i_aspec}combined', ASPEC_COLORS[i_aspec], visible=vis)
                temp_trace_count += tc

            # for d in np.unique(detectors):
            #     fig = add_tts_regions(fig, d)

            ## plot the preview spec over both of them
            fig, ntraces = plot_preview(fig, targ_dir, visible=vis)
            temp_trace_count += ntraces # preview added to each row

            n_all_traces.append(temp_trace_count) # counting the total number of traces for the buttons


        ## formatting
        fig = update_plot_layouts(fig, grating_ymaxes, target)

        fig = make_buttons(fig, all_vis, n_all_traces, buttons)

        ## save the plot out
        dr_outdir = os.path.join(outdir_name, f"dr{dr}")
        if not os.path.exists(dr_outdir):
            os.mkdir(dr_outdir)
            os.mkdir(os.path.join(dr_outdir, 'highmass'))
            os.mkdir(os.path.join(dr_outdir, 'lowmass'))

        if target in list(highmass_df['target_name_hlsp']):
            fig_outname = os.path.join(dr_outdir, f'highmass/{target}_preview_dr7.html')
        elif target in list(lowmass_df['target_name_hlsp']):
            fig_outname = os.path.join(dr_outdir, f'lowmass/{target}_preview_dr7.html')
        else:
            print(f'{target} not recognized')
            fig_outname = os.path.join(dr_outdir, f'{target}_preview_dr7.html')

        fig.write_html(fig_outname)
        print(f'Saved: {fig_outname}')

    return no_previews

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    # Trying to make more color-blind accesible. Splitting into blue & red hues
    CSPEC_COLORS = ["#125A56", "#238F9D", "#60BCE9", "#05457D", #"#C6DBED",
                    "#00767B", "#42A7C6", "#9DCCEF", "#DEE6E7"] # blues
    ASPEC_COLORS = ["#A01813", "#E94C1F", "#FD9A44", "#F9D576",
                    "#D11807", "#F57634", "#FFB954", "#F0E6B2"] # reds

    outdir = '/astro/ullyses/preview_plots/' # sorted into dr#/<highmass/lowmass>
    dr = 7 # data release number
    no_previews = make_plots(outdir, dr)
    print('No files for:')
    for p in no_previews:
        print(f"{p}; {glob.glob(os.path.join(p, '*'))}")
