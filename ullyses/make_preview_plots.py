import argparse
import glob
import os
import numpy as np
import plotly.graph_objects as go
import astropy.units as u
from plotly.subplots import make_subplots
from itertools import repeat

from astropy.io import fits
from copy import deepcopy
from specutils import Spectrum1D, SpectralRegion
from specutils.manipulation import extract_region

from ullyses_utils.match_aliases import match_aliases
from ullyses_utils.parse_csv import parse_database_csv

#-------------------------------------------------------------------------------

def make_buttons(fig, button_names):
    '''Create buttons for plotly figures, being careful to count all of the
    number of traces being added.

    Inputs:
    -------
    - fig : plotly figure
        the figure with all traces added
    - button_names : list
        What should be read on the buttons

    Outputs:
    --------
    - fig : plotly figure
        the figure returned with updated buttons
    '''

    buttons = [dict(label=str(label),
                    method='relayout',
                    args=[{'yaxis.type' : scale_type,
                           'yaxis2.type' : scale_type,
                           },])
               for label, scale_type in zip(button_names, ['linear', 'log'])]

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
        flux = flux[20000:-1001]
    elif len(flux) > 3000:
        flux = flux[1000:-1001]
    else:
        flux = flux[100:-100]

    return flux.max()

#-------------------------------------------------------------------------------

def open_files_and_plot(fig, cspec_file, row, legendgroup, color):

    #print(cspec_file)

    with fits.open(cspec_file) as hdu:
        detector = hdu[0].header['DETECTOR']
        all_gratings = [g for g in hdu[2].data['DISPERSER']]
        lbl = f"{hdu[0].header['TELESCOP']} {hdu[0].header['INSTRUME']} {'&'.join(np.unique(all_gratings))}"

        for ext in range(len(hdu[1].data)):
            wave = hdu[1].data["WAVELENGTH"][ext]
            flux = hdu[1].data["FLUX"][ext]
            # set the zeros to nans for better plotting
            flux[np.where(flux == 0.0)[0]] = "nan"

            if ext > 1:
                showlegend = False
            else:
                showlegend = True

            ## add the spectrum
            fig.add_trace(go.Scattergl(x=wave,
                                       y=flux,
                                       mode='lines',
                                       line=dict(color=color,
                                                 width=2),
                                       opacity=0.5,
                                       showlegend=False,
                                       legendgroup=legendgroup,
                                       hovertemplate = '%{text}<br>(%{x}, %{y})<extra></extra>',
                                       text = [f"<b>{'/'.join(lbl.split(' '))}</b>"]*len(wave),
                                       ), row=row, col=1)

    ## adding a trace for the label to make it bigger
    fig.add_trace(go.Scattergl(x=[wave.min()], y=[-1E-16],
                               mode='lines',
                               line=dict(color=color, width=4),
                               opacity=0.8,
                               name=lbl,
                               legendgroup=legendgroup,
                               ), row=row, col=1)

    return fig, detector

#-------------------------------------------------------------------------------

def plot_preview(fig, targ_dir):

    prev = glob.glob(os.path.join(targ_dir, '*_preview-spec.fits'))
    if len(prev) == 0:
        #setting the default no preview to something reasonable
        flux_max = [5E-13]
        return fig, flux_max

    with fits.open(prev[0]) as hdu:
        wave = hdu[1].data["WAVELENGTH"].ravel()
        flux = hdu[1].data["FLUX"].ravel()

        echelle = any([g.startswith('E') for g in hdu[2].data['DISPERSER']])
        # adding in a deepcopy to the flux b/c I modify the array & python edits it and
        #   then doesn't plot it otherwise
        flux_max = get_max_flux(wave, deepcopy(flux), echelle=echelle)

        # after getting the plotting maxes, set the zero flux to nan for better plots
        flux[np.where(flux == 0.0)[0]] = "nan"

    for row in [1, 2]:
        fig.add_trace(go.Scatter(x=wave,
                                 y=flux,
                                 mode='lines',
                                 line=dict(width=2,
                                           color='dimgrey'),
                                 opacity=0.6,
                                 legendgroup='preview',
                                 showlegend=False,
                                 hovertemplate = '<b>All Abutted</b><br>(%{x}, %{y})<extra></extra>',
                                 legendrank=1001, # in legend top if legend rank < 1000 & bottom if rank > 1000
                                 ), row=row, col=1)

    ## adding a trace for the label to make it bigger
    fig.add_trace(go.Scattergl(x=[wave.min()], y=[-1E-16],
                               mode='lines',
                               line=dict(color='dimgrey', width=4),
                               opacity=0.8,
                               name='All Abutted Spectra<br>(Level 4 HLSP)',
                               legendgroup='preview',
                               ), row=1, col=1)

    return fig, flux_max

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
                             'range' : [-0.1*np.max(ymax_arr), 1.10*np.max(ymax_arr)],
                             'title': {'text' : 'Flux (erg/s/cm<sup>2</sup>/&#8491;)',
                                       'font' : {'size' : 25}, },
                             'tickfont' : {'size' : 20},
                             },
                      yaxis2={'showexponent': 'all',
                              'exponentformat': 'E',
                              'autorange' : False,
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
                      legend={'font' : {'size' : 15} },
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

def make_plots(base_datadir, outdir_name, dr):

    # There are different levels for different combinations of gratings:
    #   panel 1 is individual gratings
    #   panel 2 is combined gratings
    #   overplotted on both is the "preview_spec" which can be turned on&off

    highmass_df, lowmass_df = high_and_low_df()

    # v-cv-cha
    # hd-104237e
    # av-456
    no_previews = []
    for targ_dir in np.sort(glob.glob(os.path.join(base_datadir, "*"))):

        targ_dir = os.path.join(targ_dir, f"dr{dr}")

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

        ## row 1; grating only files (cspec)
        detectors = []
        for i_cspec, cspec_file in enumerate(np.sort(grating_cspec)):
            fig, d = open_files_and_plot(fig, cspec_file, 1, f'{i_cspec}grating', CSPEC_COLORS[i_cspec])
            detectors.append(d)

        ## row 2; abutted files
        for i_aspec, aspec_file in enumerate(np.sort(combined_aspec)):
            fig, d = open_files_and_plot(fig, aspec_file, 2, f'{i_aspec}combined', ASPEC_COLORS[i_aspec])

        # for d in np.unique(detectors):
        #     fig = add_tts_regions(fig, d)

        ## plot the preview spec over both of them
        fig, plotting_max = plot_preview(fig, targ_dir)

        ## formatting
        fig = update_plot_layouts(fig, plotting_max, target)

        # create empty buttons for multiple scale options
        buttons = ['Linear Scale', 'Log Scale']
        fig = make_buttons(fig, buttons)

        ## save the plot out
        dr_outdir = os.path.join(outdir_name, f"dr{dr}")
        if not os.path.exists(dr_outdir):
            os.mkdir(dr_outdir)
            os.mkdir(os.path.join(dr_outdir, 'highmass'))
            os.mkdir(os.path.join(dr_outdir, 'lowmass'))

        if target in list(highmass_df['target_name_hlsp']):
            fig_outname = os.path.join(dr_outdir, f'highmass/{target}_preview_dr{dr}.html')
        elif target in list(lowmass_df['target_name_hlsp']):
            fig_outname = os.path.join(dr_outdir, f'lowmass/{target}_preview_dr{dr}.html')
        else:
            print(f'{target} not recognized')
            fig_outname = os.path.join(dr_outdir, f'{target}_preview_dr{dr}.html')

        fig.write_html(fig_outname)
        print(f'Saved: {fig_outname}')

    return no_previews

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    # Trying to make more color-blind accesible. Splitting into blue & red hues
    # ASPEC_COLORS = ["#125A56", "#238F9D", "#60BCE9", "#05457D", #"#C6DBED",
    #                 "#00767B", "#42A7C6", "#9DCCEF", "#DEE6E7"] # blues
    ASPEC_COLORS = ["#C2A5CF", "#9970AB", "#762A83", "#E7D4E8"] # purple
    CSPEC_COLORS = ["#A01813", "#D11807", "#E94C1F", "#F57634", "#FD9A44",
                    "#FFB954"] # reds
    # CSPEC_COLORS = CSPEC_COLORS[::-1]
    # ASPEC_COLORS = ASPEC_COLORS[::-1]

    datadir = "/astro/ullyses/ULLYSES_HLSP/"
    outdir = '/astro/ullyses/preview_plots/' # sorted into dr#/<highmass/lowmass>
    dr = 7 # data release number
    no_previews = make_plots(datadir, outdir, dr)
    print('No files for:')
    for p in no_previews:
        print(f"{p}; {glob.glob(os.path.join(p, '*'))}")
