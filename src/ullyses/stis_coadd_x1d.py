import shutil
import argparse
import glob
import os
import numpy as np
from astropy.io import fits

from ullyses.ullyses_coadd_abut_wrapper import Ullyses_STISSegmentList
from ullyses.fuse_add_dq import add_column

"""
This code takes a list of STIS x1ds and coadds
the spectra over a uniform wavelength grid.

Arguments:
    files (list): Space-separated list of files to coadd
    targ (str): Target name
    outdir (str): Path to output directory
"""


class STIScoadd(Ullyses_STISSegmentList):
    """
    A class to perform coadditions of STIS x1ds.
    """

    def create_output_wavelength_grid(self):
        """
        Using the input spectra, create a common wavelength grid
        over which to perform the coaddition.
        :return: wavelength grid
        """

        min_wavelength = 10000.0
        max_wavelength = 0.0
        for segment in self.members:
            minwave = segment.data['wavelength'].min()
            maxwave = segment.data['wavelength'].max()
            if minwave < min_wavelength: min_wavelength = minwave
            if maxwave > max_wavelength: max_wavelength = maxwave

        max_delta_wavelength = 0.0
        for segment in self.members:
            wavediffs = segment.data['wavelength'][1:] - segment.data['wavelength'][:-1]
            max_delta_wavelength = max(max_delta_wavelength, wavediffs.max())

        self.delta_wavelength = max_delta_wavelength

        self.min_wavelength = int(min_wavelength)
        self.max_wavelength = int(max_wavelength+self.delta_wavelength) + 1

        wavegrid = np.arange(self.min_wavelength, self.max_wavelength, self.delta_wavelength)

        self.output_wavelength = wavegrid
        self.nelements = len(wavegrid)
        self.output_sumflux = np.zeros(self.nelements)
        self.output_sumweight = np.zeros(self.nelements)
        self.output_flux = np.zeros(self.nelements)
        self.output_errors = np.zeros(self.nelements)
        self.signal_to_noise = np.zeros(self.nelements)
        self.output_exptime = np.zeros(self.nelements)

        return wavegrid

    def coadd(self, ignore_dq_file):
        """
        Coadd the spectra.
        :param ignore_dq_file: Fits file with DQ array to remove from coadd
        :return: None
        """

        self.output_dq = np.zeros(self.nelements).astype(int)
        self.output_gross = np.zeros(self.nelements)
        self.output_net = np.zeros(self.nelements)
        self.output_sumgross = np.zeros(self.nelements)
        self.output_sumnet = np.zeros(self.nelements)
        self.output_sumerrors = np.zeros(self.nelements)
        for i, segment in enumerate(self.members):
            if os.path.basename(self.datasets[i]) == ignore_dq_file:
                goodpixels = np.arange(len(segment.data['dq']))
            else:
                goodpixels = np.where((segment.data['dq'] & segment.sdqflags) == 0)
            wavelength = segment.data['wavelength'][goodpixels]
            indices = self.wavelength_to_index(wavelength)
            gross_counts = self.get_flux_weight(segment)
            weight = gross_counts[goodpixels]
            flux = segment.data['flux'][goodpixels]
            gross = segment.data['gross'][goodpixels]
            net = segment.data['net'][goodpixels]
            err = segment.data['error'][goodpixels]
            self.output_sumerrors[indices] = np.sqrt((self.output_sumerrors[indices]**2) + (err**2))
            self.output_sumgross[indices] += gross
            self.output_sumnet[indices] += net
            if os.path.basename(self.datasets[i]) != ignore_dq_file:
                all_indices = self.wavelength_to_index(segment.data['wavelength'])
                self.output_dq[all_indices] = self.output_dq[all_indices] | segment.data['dq']
            self.output_sumweight[indices] = self.output_sumweight[indices] + weight
            self.output_sumflux[indices] = self.output_sumflux[indices] + flux * weight
            self.output_exptime[indices] = self.output_exptime[indices] + segment.exptime                         
        good_dq = np.where(self.output_exptime > 0.)
        self.first_good_wavelength = self.output_wavelength[good_dq][0]
        self.last_good_wavelength = self.output_wavelength[good_dq][-1]
        nummembers = len(self.members)
        nonzeros = np.where(self.output_sumweight == nummembers)
        self.output_flux[nonzeros] = self.output_sumflux[nonzeros]
        self.output_gross[nonzeros] = self.output_sumgross[nonzeros]
        self.output_net[nonzeros] = self.output_sumnet[nonzeros]
        self.output_errors[nonzeros] = self.output_sumerrors[nonzeros]
        self.signal_to_noise[nonzeros] = self.output_sumweight[nonzeros] / self.output_errors[nonzeros]
        return


def coadd_1d_spectra(files, targ, outdir):
    """
    Perform the coaddition of the x1ds and
    format the header of the output fits file.
    :param files: Space-separated list of files to coadd
    :param targ: Target name
    :param outdir: Path to output directory
    :return: None
    """

    targ = targ.upper()
    grating = fits.getval(files[0], "opt_elem")
    root = fits.getval(files[0], "rootname").lower()
    coadd_dir = os.path.join(outdir, f"{grating}_coadd")
    if not os.path.exists(coadd_dir):
        os.makedirs(coadd_dir)
    for item in files:
        shutil.copy(item, coadd_dir)
        print(f"Copied {item} to {coadd_dir}")
    combined0 = f"{root}_{targ.lower()}_x1d.fits"
    combined = os.path.join(outdir, combined0)
    if os.path.exists(combined):
        os.remove(combined)
        print(f"Removed {combined}")
    prod = STIScoadd(grating, inpath=coadd_dir, weighting_method='unity')
    prod.target = prod.get_targname()
    prod.targ_ra, prod.targ_dec, prod.coord_epoch = prod.get_coords()
    prod.create_output_wavelength_grid()
    ignore_file = os.path.basename(files[1])
    prod.coadd(ignore_dq_file=ignore_file)

    nelements = len(prod.output_net)
    wl_arr = prod.output_wavelength.reshape(1, nelements)
    flux_arr = prod.output_flux.reshape(1, nelements)
    err_arr = prod.output_errors.reshape(1, nelements)
    net_arr = prod.output_net.reshape(1, nelements)
    gross_arr = prod.output_gross.reshape(1, nelements)
    dq_arr = prod.output_dq.reshape(1, nelements)
    hdr0 = fits.getheader(files[0], 0)
    hdr0["TARGNAME"] = targ
    new_hdu0 = fits.PrimaryHDU(header=hdr0)
    cols = []
    cols.append(fits.Column(name="WAVELENGTH", format=f"{nelements}D",
                            array=wl_arr, unit="Angstroms"))
    cols.append(fits.Column(name="FLUX", format=f"{nelements}E",
                            array=flux_arr, unit="erg/s/cm**2/Angstrom"))
    cols.append(fits.Column(name="ERROR", format=f"{nelements}E",
                            array=err_arr, unit="erg/s/cm**2/Angstrom"))
    cols.append(fits.Column(name="NET", format=f"{nelements}E",
                            array=net_arr, unit="Counts/s"))
    cols.append(fits.Column(name="GROSS", format=f"{nelements}E",
                            array=gross_arr, unit="Counts/s"))
    cols.append(fits.Column(name="DQ", format=f"{nelements}I",
                            array=dq_arr, unit=None))
    coldefs = fits.ColDefs(cols)                                                                            
    hdr = fits.getheader(files[0], 1)                                                                   
    t = fits.BinTableHDU.from_columns(coldefs, header=hdr)                                                  
    final_hdus = [new_hdu0, t]                                                                            
    new_hdulist = fits.HDUList(final_hdus)                                                                  
    new_hdulist.writeto(combined, overwrite=True)                                                         
    print(f"Wrote {combined}")                                                                            
    # add_column(combined, combined, 1, "WAVELENGTH", f"{nelements}D", wl_arr, colunit="Angstroms", overwrite=True)
    # add_column(combined, combined, 1, "FLUX", f"{nelements}E", flux_arr, colunit="erg/s/cm**2/Angstrom", overwrite=True)
    # add_column(combined, combined, 1, "ERROR", f"{nelements}E", err_arr, colunit="erg/s/cm**2/Angstrom", overwrite=True)
    # add_column(combined, combined, 1, "NET", f"{nelements}E", net_arr, colunit="Counts/s", overwrite=True)
    # add_column(combined, combined, 1, "GROSS", f"{nelements}E", gross_arr, colunit="Counts/s", overwrite=True)
    # add_column(combined, combined, 1, "DQ", f"{nelements}I", dq_arr, overwrite=True)


def parse_input(files):
    """
    Ensure the input files are in the correct format.
    :param files: Space-separated list of files to coadd
    :return: None
    """

    if len(files) == 1 or "," in files[0]:
        raise TypeError(f"Names of input files must be separated by spaces. Do not include commas.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", nargs="+",
                        help="Space-separated list of files to coadd")
    parser.add_argument("-t", "--targ", 
                        help="Name of output target name")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Output directory")
    args = parser.parse_args()
    parse_input(args.files)
    coadd_1d_spectra(args.files, args.targ, args.outdir) 
