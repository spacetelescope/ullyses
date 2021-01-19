#! /usr/bin/env python

import datetime
import argparse
import os
from astropy.io import fits
import shutil
import glob
import itertools
import numpy as np
import sys
import stistools

from read_config import read_config
os.environ["oref"] = "/grp/hst/cdbs/oref/"

class Stisdata():
    """
    A class to describe STIS data throughout the calibration process.

    Attributes:
        outdir (str): Path where output products should be written.
        infiles (list): List of input raw datasets.
        nfiles (int): Number of input raw datasets.
        basedir (str): Path to raw STIS data.
        rootname (dict): Dictionary where each key is a dictionary- each key
            is the rootname and vals describe dataset properties (detector,
            dither offset, actual file name).
        extract_info (dict): Nested dictionary where keys are rootnames; each
            nested dictionary describes dataset properties (visit, detector,
            XD shift, raw filename, flt filename, shifted flt filename, and
            combined image filename.
        _sci_dir (str): 'sci' directory used by stis_cti
        _dark_dir (str): 'dark' directory used by stis_cti
        _ref_dir (str): 'ref' directory used by stis_cti
#!!!        combined (dict): Nested dictionary where keys are the filenames of 
            the combined FLTs and FLCs; each nested dict listing the
            visit and detector values.
        x1d (list): List of final x1d products. 
        
    """
    def __init__(self, scifile, outdir=None, yamlfile=None):
        """
        Args:
            infiles: List or wildcard describing all input files.
            outdir: Output directory for products.
        """

        detector = fits.getval(scifile, "detector")
        opt_elem = fits.getval(scifile, "opt_elem")
        assert detector == "CCD" and opt_elem !="MIRVIS", f"Observing configurtion not supported: {detector}/{opt_elem}"
        self.scifile = scifile
        self.basedir = os.path.dirname(self.scifile)
        if outdir is None:
            nowdt = datetime.datetime.now()
            outdir = os.path.join(self.basedir, nowdt.strftime("%Y%m%d_%H%M"))
        self.outdir = outdir
        self.rootname = fits.getval(scifile, "rootname")
        self.visit = self.rootname[4:6]
        self.target = fits.getval(scifile, "targname")
        if yamlfile is None:
            self.x1d_c, self.crrej_c, self.defringe_c, self.cti_proc = read_config(target=self.target)
        else:
            self.x1d_c, self.crrej_c, self.defringe_c, self.cti_proc = read_config(yamlfile=yamlfile)
        self.fringeflat = self.defringe_c["fringeflat"]
        self.opt_elem = opt_elem


#-----------------------------------------------------------------------------#

    def run_all(self):
        self.perform_cti()
        self.analyze_dark()
        self.flag_negatives()
        self.crrej()
        self.defringe()
        self.extract_spectra()

#-----------------------------------------------------------------------------#

    def flag_negatives(self, change_val=False, dq=4, thresh=-100):
        """
        For pixels with values below the input negative threshold, change the DQ
        value to a specified value. Also, if specified, change large negative 
        values to large positive values.
    
        Args:
            change_val (Bool): True if large negative pixel values should
                be chnaged to large positive values, default is False.
            dq (int): DQ value that large negative pixel values should
                be assigned, default is 4.
            thresh (int or float): Threshold below wich pixels should
                be flagged with specified DQ flag, default is -100.
        Returns:
            None
        """

        print(f"Flagging pixels with counts below {thresh} with DQ={dq}...")
        with fits.open(self.flt, mode="update") as sci_hdu: 
            sci_data = sci_hdu[1].data
            neg = np.where(sci_data <= thresh)
            sci_hdu[3].data[neg] += dq
            if change_val:
                print("Changing large negative values to large positive values...")
                sci_hdu[1].data[neg] = 10000

#-----------------------------------------------------------------------------#

    def analyze_dark(self, dq=16):
        """
        Create custom superdarks that are identical to those in FLT headers'
        except that they have correct dark-flagging done (STIS reference 
        pipeline does not do it correctly). Apply dark DQ (by default, DQ=16: 
        "Pixel having dark rate >5σ times the median dark level" -STIS IHB) 
        to FLTs.

        Args:
            dq (int): DQ to apply to pixels with large dark values.
        Returns:
            None 
        """

        print(f"Manually creating superdarks and setting DQ={dq} values...")
        customdark_dir = os.path.join(self.basedir, "custom_darks")
        if not os.path.isdir(customdark_dir):
            os.mkdir(customdark_dir)
            print(f"Made directory: {customdark_dir}")
            
        # Read in science FLT dataset.
        sci_hdu = fits.open(self.flt, mode="update")
        sci_dq = sci_hdu[3].data
        darkfile0 = sci_hdu[0].header["darkfile"]
    
        # Determine DARKFILE filename.
        if "/" in darkfile0:
            darkname = darkfile0.split("/")[1]
            darkfile = os.path.join(self._ref_dir, darkname)
        else:
            darkname = darkfile0.split("$")[1]
            darkfile = os.path.join(OREF_DIR, darkname)

        # Remove any existing DQ={dq} flags from FLT, since DQ={dq} flags are
        # inherently wrong.
        outfile = os.path.join(customdark_dir, darkname)
        sci_hdu[0].header["DARKFILE"] = outfile
        # DQ=512 is "bad pixel in reference file"
        sci_dq16 = np.where((sci_dq&dq == dq) & (sci_dq&512 == 0))
        sci_hdu[3].data[sci_dq16] -= dq
        
        written_darks = [os.path.basename(x) for x in 
                         glob.glob(os.path.join(customdark_dir, "*fits"))]
        # If the custom darkfile has already been written (from a previous iteration), 
        # read from the darkfile and apply the DQ={dq} flags to the FLT. 
        if darkname in written_darks:
            dark_dq = fits.getdata(outfile, 3)
            dark_dq16 = np.where(dark_dq&dq == dq)
            sci_hdu[3].data[dark_dq16] |= dq
            sci_hdu.close()
            print(f"Already wrote {darkname}")
        
        # Create custom darkfile, flag high dark values with dq={dq} 
        else:
            dark_hdu = fits.open(darkfile)
            dark = dark_hdu[1].data
            dark_dq = dark_hdu[3].data
            dark_dq16 = np.where(dark_dq&dq == dq)
            dark_dq[dark_dq16] -= dq
            sci_hdu[3].data[dark_dq16] |= dq
        
            # Determine 5*stddev + median of the darkfile, anything above this is DQ={dq}. 
            dark_thresh = (np.std(dark)*5) + np.median(dark)
            dark_inds = np.where((dark > dark_thresh) | (dark < -dark_thresh))
            dark_dq[dark_inds] |= dq
        
            new_dark = fits.HDUList( fits.PrimaryHDU() )
            new_dark[0].header = dark_hdu[0].header
            for ext in [1,2]:
                new_dark.append(fits.ImageHDU(dark_hdu[ext].data, dark_hdu[ext].header, name=dark_hdu[ext].name))
            new_dark.append(fits.ImageHDU(dark_dq, dark_hdu[3].header, name=dark_hdu[ext].name))
            new_dark.writeto(outfile)
        
            sci_hdu.close()
            dark_hdu.close()
        
            print(f"Wrote new darkfile: {outfile}")
                    
#-----------------------------------------------------------------------------#

    def perform_cti(self):
        """
        Run the STIS CTI code on STIS CCD data.
        """
        
        import stis_cti
        print("Performing CTI correction on data...") 
        
        self._sci_dir = os.path.join(self.basedir, "science")
        self._dark_dir = os.path.join(self.basedir, "darks")
        self._ref_dir = os.path.join(self.basedir, "ref")
        
        # These directories need to exist for stis_cti to run. 
        for direc in [self._dark_dir, self._ref_dir, self._sci_dir]:
            if not os.path.exists(direc):
                os.mkdir(direc)
                print(f"Made directory: {direc}")

        # stis_cti needs raw, epc, spt, asn, and wav files as input.
        print("Copying CCD datasets to science directory...")
        allfiles = glob.glob(os.path.join(self.basedir, "*_raw.fits")) + \
                   glob.glob(os.path.join(self.basedir, "*_epc.fits")) + \
                   glob.glob(os.path.join(self.basedir, "*_spt.fits")) + \
                   glob.glob(os.path.join(self.basedir, "*_asn.fits")) + \
                   glob.glob(os.path.join(self.basedir, "*_wav.fits"))
        for item in allfiles:
            shutil.copy(item, self._sci_dir)
    
        # Run stis_cti
        stis_cti.stis_cti(self._sci_dir, self._dark_dir, self._ref_dir, self.cti_proc, verbose=True)

        self.copy_products()
        
#-----------------------------------------------------------------------------#

    def copy_products(self):
        """
        Copy the intermediate products from STIS CTI products to self.outdir, 
        since several calibration iterations are often needed and we don't want
        to overwrite the original products.
        """

        if os.path.exists(self.outdir):
            print("WARNING: Output directory already exists, deleting {0}".format(self.outdir))
            shutil.rmtree(self.outdir)
        os.mkdir(self.outdir)
    
        # If stis_cti was not run, define the stis_cti products directories.
        if not hasattr(self, "_sci_dir"):
            self._sci_dir = os.path.join(self.basedir, "science")
            self._dark_dir = os.path.join(self.basedir, "darks")
            self._ref_dir = os.path.join(self.basedir, "ref")
       
        # For CCD data, copy FLCs. For MAMA data, copy FLTs.
        flc_files = glob.glob(os.path.join(self._sci_dir, "*flc.fits"))
        for item in flc_files:
            shutil.copy(item, self.outdir)
        
        print(f"Copied FLC files to {self.outdir}")
        self.flt = os.path.join(self.outdir, self.rootname+"_flc.fits")
        
#-----------------------------------------------------------------------------#

    def extract_spectra(self):
        """
        Extract the spectra using stistools.x1d.
        """

        from stistools import x1d
        print("Extracting spectra...")
        
        # Look up extraction parameters given the target config file.
        for targ, pars in self.x1d_c.items():
            if targ == "sci":
                outfile = os.path.join(self.outdir, self.rootname+"_x1d.fits")
                sci_x1d = outfile
            else:
                outfile = os.path.join(self.outdir, f"{self.rootname}_{targ}_x1d.fits")
            x1d.x1d(self.drj, 
                output = outfile, 
                a2center = pars["yloc"], 
                maxsrch = pars["maxsrch"],
                extrsize = pars["height"],
                bk1offst = pars["b_bkg1"] - pars["yloc"],
                bk2offst = pars["b_bkg2"] - pars["yloc"],
                bk1size = pars["b_hgt1"],
                bk2size = pars["b_hgt2"],
                ctecorr="omit",
                verbose=True)
            print(f"Wrote x1d file: {outfile}")

        self.x1d = sci_x1d

#-----------------------------------------------------------------------------#

    def crrej(self):
        from stistools.ocrreject import ocrreject
        print("Performing cosmic ray rejection...")
        
        outfile = os.path.join(self.outdir, self.rootname+"_crc.fits")
        
        # Look up crrej parameters given the target config file.
        ocrreject(self.flt,
                  output=outfile,
                  initgues=self.crrej_c["initgues"],
                  crsigmas=self.crrej_c["crsigmas"],
                  crradius=self.crrej_c["crradius"],
                  crthresh=self.crrej_c["crthresh"],
                  crmask=self.crrej_c["crmask"],
                  verbose=self.crrej_c["verbose"])

        print(f"Wrote crj file: {outfile}")
        self.crj = outfile

#-----------------------------------------------------------------------------#

    def defringe(self):

        fringeflat = self.defringe_c["fringeflat"]
        rawfringe = os.path.join(self.basedir, fringeflat)
        fringeroot = os.path.basename(fringeflat)[:9]
        outnorm = os.path.join(self.outdir, fringeroot+"_nsp.fits")
        stistools.defringe.normspflat(inflat=rawfringe,
                   do_cal=True,
                   outflat=outnorm,
                   biasfile=self.defringe_c["normspflat"]["biasfile"],
                   darkfile=self.defringe_c["normspflat"]["darkfile"])
        outmk = os.path.join(self.outdir, fringeroot+"_mff.fits")
        stistools.defringe.mkfringeflat(inspec=self.crj,
                     inflat=outnorm,
                     outflat=outmk,
                     do_shift=self.defringe_c["mkfringeflat"]["do_shift"],
                     beg_shift=self.defringe_c["mkfringeflat"]["beg_shift"],
                     end_shift=self.defringe_c["mkfringeflat"]["end_shift"],
                     shift_step=self.defringe_c["mkfringeflat"]["shift_step"],
                     do_scale=self.defringe_c["mkfringeflat"]["do_scale"],
                     beg_scale=self.defringe_c["mkfringeflat"]["beg_scale"],
                     end_scale=self.defringe_c["mkfringeflat"]["end_scale"],
                     scale_step=self.defringe_c["mkfringeflat"]["scale_step"])
        outfile = stistools.defringe.defringe(science_file=self.crj,
                           fringe_flat=outmk,
                           verbose=True)

        print(f"Wrote defringed crj file: {outfile}")
        self.drj = outfile
        
