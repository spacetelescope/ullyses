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

def add_spectral_trace(fig, wave, flux, lbl, row, legendgroup, showlegend=True, dash=None, visible=True):

    fig.add_trace(go.Scattergl(x=wave,
                               y=flux,
                               mode='lines',
                               line=dict(dash=dash,
                                         color=COLORS[i],
                                         width=4),
                               opacity=0.8,
                               name=lbl,
                               showlegend=showlegend,
                               legendgroup=legendgroup,
                               visible=visible,
                               hovertemplate = '%{text}<br>(%{x}, %{y})<extra></extra>',
                               text = [f"<b>{'/'.join(lbl.split(' '))}</b>"]*len(wave),
                               ), row=row, col=1)

    return fig

#-------------------------------------------------------------------------------

def open_files_and_plot(fig, cspec_file, row, flux_maxes, legendgroup, visible=True, dash=None):

    print(cspec_file)
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

            fig = add_spectral_trace(fig, wave, flux, lbl, row, legendgroup,
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

    for row in [1, 2]:

        if row == 1:
            showlegend = True
        else:
            showlegend = False

        fig.add_trace(go.Scatter(x=wave,
                                 y=flux,
                                 mode='lines',
                                 line=dict(width=4, color='dimgrey'),
                                 opacity=0.7,
                                 name='All Abutted Spectra<br>(Level 4 HLSP)',
                                 legendgroup='preview',
                                 showlegend=showlegend,
                                 visible=visible,
                                 hovertemplate = '<b>All Abutted</b><br>(%{x}, %{y})<extra></extra>',
                                 ), row=row, col=1)

        ntraces += 1

    return fig, ntraces

#-------------------------------------------------------------------------------

def update_plot_layouts(fig, ymax_arr):

    fig.update_annotations(font_size=25)

    fig.add_annotation(text="<sup>*</sup>Gratings with the same resolution ",
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

if __name__ == "__main__":

    # open all of the non-combined(?) cspec files for a specific target
    #   this will be *_<target>_<grating>_dr?_cspec.fits
    #   compared to the grating combined ones that will have "-s" in them
    # This is what Alex does in his plotting routine

    # Maybe we could have different panels for the different levels of combinations?
    #   panel 1 is individual gratings
    #   panel 2 is combined gratings
    #   overplotted on both is the "preview_spec" which can be turned on&off
    #     -> or this can be in panel 3

    COLORS = ["#9970ab", "#5aae61", "#d95f02", "#e7298a", "#1bddf2", "#1f78b4", "#fb9a99",
              "#fdbf6f", "#ffe44b", "#b15928", "#cab2d6", "#b2df8a", "#749688", "#7a7a7a",
              "#911cff", "#a6cee3"]

    # v-cv-cha
    # hd-104237e
    # av-456
    targ_dir = "/astro/ullyses/ULLYSES_HLSP/av-456/dr6/"
    all_cspec_files = glob.glob(os.path.join(targ_dir, '*_cspec.fits'))

    ## sort out the different types of cspec files
    grating_cspec = []
    combined_cspec = []
    for f in all_cspec_files:
        temp_grating = fits.getval(f, 'HLSP_LVL')
        if temp_grating == 3:
            # grating combined spectra
            combined_cspec.append(f)
        else:
            # potentially cenwave combined, but only 1 grating
            grating_cspec.append(f)

    # initialize a file

    #fig = go.Figure()
    fig = make_subplots(rows=2,
                        cols=1,
                        shared_xaxes=True,
                        #shared_yaxes=True, # this only shares across columns...
                        subplot_titles=('Single Grating Coaddition',
                                        'Single Instrument, Multiple Grating<sup>*</sup> Abutment'),
                        row_titles=('<a href="https://ullyses.stsci.edu/ullyses-data-description.html#productDescrip">Level 2 HLSP</a>',
                                    '<a href="https://ullyses.stsci.edu/ullyses-data-description.html#productDescrip">Level 3 HLSP</a>'),
                        vertical_spacing=0.05)

    target = fits.getval(grating_cspec[0], 'TARGNAME')

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

        ## row 1; grating only files
        grating_ymaxes = []
        detectors = []
        for i, cspec_file in enumerate(np.sort(grating_cspec)):
            fig, grating_ymaxes, d, tc = open_files_and_plot(fig, cspec_file, 1, grating_ymaxes, f'{i}grating', visible=vis)
            detectors.append(d)
            temp_trace_count += tc

        ## row 2; abutted files
        combined_maxes = []
        for i, cspec_file in enumerate(np.sort(combined_cspec)):
            fig, combined_maxes, d, tc = open_files_and_plot(fig, cspec_file, 2, combined_maxes, f'{i}combined', dash='dot', visible=vis)
            temp_trace_count += tc

        # for d in np.unique(detectors):
        #     fig = add_tts_regions(fig, d)

        ## plot the preview spec over both of them
        fig, ntraces = plot_preview(fig, targ_dir, visible=vis)
        temp_trace_count += ntraces # preview added to each row

        n_all_traces.append(temp_trace_count) # counting the total number of traces for the buttons


    ## formatting
    fig = update_plot_layouts(fig, grating_ymaxes)

    fig = make_buttons(fig, all_vis, n_all_traces, buttons)


    fig.show()
    fig.write_html(f'/Users/rplesha/Desktop/{target}_preview_v4.html')
