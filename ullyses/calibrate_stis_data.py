#! /usr/bin/env python
from __future__ import print_function, division

"""
Unfortunately, ensure you are in the IRAF27 environment before running this scripts.
"""

import argparse
import os
from astropy.io import fits
import shutil
import glob
import itertools
import numpy as np
import sys

from parameters import *

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

class stisdata():
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
        combined (dict): Nested dictionary where keys are the filenames of 
            the combined FLTs and FLCs; each nested dict listing the
            visit and detector values.
        x1d (list): List of final x1d products. 
        
    """
    def __init__(self, infiles, outdir):
        """
        Args:
            infiles: List or wildcard describing all input files.
            outdir: Output directory for products.
        """

        self.outdir = outdir
        self.infiles = self.parse_infiles(infiles)
        self.nfiles = len(self.infiles)
        self.basedir = os.path.dirname(self.infiles[0])
        self.rootname = {}

        # Exclude MIRVIS, which are not science data.
        for rawfile in self.infiles:
            if fits.getval(rawfile, "opt_elem") != "MIRVIS":
                rootname = fits.getval(rawfile, "rootname")
                self.rootname[rootname] = {}
                self.rootname[rootname]["visit"] = rootname[4:6]
                detector = fits.getval(rawfile, "detector")
                if "MAMA" in detector:
                    self.rootname[rootname]["detector"] = "mama"
                elif "CCD" in detector:
                    self.rootname[rootname]["detector"] = "ccd"
                else:
                    self.rootname[rootname]["detector"] = detector
                self.rootname[rootname]["shift"] = SHIFT_DICT[rootname]["shift"]
                self.rootname[rootname]["rawfile"] = rawfile
        
        # Imported from parameters.py
        self.extract_info = EXTRACT_INFO

#-----------------------------------------------------------------------------#

    def parse_infiles(self, infiles):
        """
        Gather all input raw files into a list.

        Args:
            infiles: List or wild-card string of input datasets.
        Returns:
            rawfiles: List of input datasets.
        """

        if isinstance(infiles, list):
            rawfiles = infiles

        elif isinstance(infiles, str):
            if "*" in infiles:
                rawfiles = glob.glob(infiles)
                assert len(rawfiles) != 0, \
                "ERROR: No matching datasets for {0}".format(infiles)
            else:
                rawfiles = infiles

        return rawfiles

#-----------------------------------------------------------------------------#

    def flag_negatives(self, change_val=False, dq=4, thresh=-100):
        """
        For pixels with values below the negative threshold, change the DQ
        value to a specified value. Also, if specified, change large negative 
        values to large positive values.
    
        Args:
            change_val (Bool): True if large negative pixel values should
                be chnaged to large positive values.
            dq (int): DQ value that large negative pixel values should
                be assigned.
            thresh (int or float): Threshold below wich pixels should
                be flagged with specified DQ flag.
        Returns:
            None
        """

        print("Flagging pixels with counts below {0} with DQ={1}...".format(
              thresh, dq))
        for rootname in self.rootname:
            flt = self.rootname[rootname]["flt"]
            with fits.open(flt, mode="update") as sci_hdu: 
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

        print("Manually creating superdarks and setting DQ={0} values...".format(dq))
        customdark_dir = os.path.join(self.basedir, "custom_darks")
        if not os.path.isdir(customdark_dir):
            os.mkdir(customdark_dir)
            print("Made directory: {0}".format(customdark_dir))

        # Do not process MAMA datasets since the dark current is so low.
        for rootname in self.rootname:
            if self.rootname[rootname]["detector"] != "ccd":
                continue
            
            # Read in science FLT dataset.
            flt = self.rootname[rootname]["flt"]
            sci_hdu = fits.open(flt, mode="update")
            sci_dq = sci_hdu[3].data
            darkfile0 = sci_hdu[0].header["darkfile"]
    
            # Determine DARKFILE filename.
            if "/" in darkfile0:
                darkname = darkfile0.split("/")[1]
                darkfile = os.path.join(self._ref_dir, darkname)
            else:
                darkname = darkfile0.split("$")[1]
                darkfile = os.path.join(OREF_DIR, darkname)

            # Remove any existing DQ=dq flags from FLT, since DQ=dq flags are
            # inherently wrong.
            outfile = os.path.join(customdark_dir, darkname)
            sci_hdu[0].header["DARKFILE"] = outfile
            sci_dq16 = np.where((sci_dq&dq == dq) & (sci_dq&512 == 0))
            sci_hdu[3].data[sci_dq16] -= dq
        
            written_darks = [os.path.basename(x) for x in 
                             glob.glob(os.path.join(customdark_dir, "*fits"))]
            # If the custom darkfile has already been written, read 
            # from the darkfile and apply the DQ=dq flags to the FLT. 
            if darkname in written_darks:
                dark_dq = fits.getdata(outfile, 3)
                dark_dq16 = np.where(dark_dq&dq == dq)
                sci_hdu[3].data[dark_dq16] |= dq
                sci_hdu.close()
                print("Already wrote {0}".format(darkname))
            
            # Create custom darkfile, flag values 
            else:
                dark_hdu = fits.open(darkfile)
                dark = dark_hdu[1].data
                dark_dq = dark_hdu[3].data
                dark_dq16 = np.where(dark_dq&dq == dq)
                dark_dq[dark_dq16] -= dq
                sci_hdu[3].data[dark_dq16] |= dq
        
                # Determine 5*stddev + median of the darkfile, anything above this is DQ=dq.    
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
            
                print("Wrote new darkfile to {0}".format(outfile))
                    
#-----------------------------------------------------------------------------#

    def perform_cti(self, processes=15):
        """
        Run the STIS CTI code on STIS CCD data.

        Args:
            processes (int): Number of processes to use while running stis_cti.
        """
        
        import stis_cti
        print("Performing CTI correction on data in {0}...".format(self.basedir))
        
        self._sci_dir = os.path.join(self.basedir, "science")
        self._dark_dir = os.path.join(self.basedir, "darks")
        self._ref_dir = os.path.join(self.basedir, "ref")
        
        # These directories need to exist for stis_cti to run. 
        for direc in [self._dark_dir, self._ref_dir, self._sci_dir]:
            if not os.path.exists(direc):
                os.mkdir(direc)
                print("Made directory: {0}".format(direc))
    
        # stis_cti needs raw, epc, spt, asn, and wav files as input.
        print("Copying CCD datasets to {0}...".format(self._sci_dir))
        for rootname in self.rootname:
            if self.rootname[rootname]["detector"] == "ccd":
                allfiles = glob.glob(os.path.join(self.basedir, rootname+"*_raw.fits")) + \
                           glob.glob(os.path.join(self.basedir, rootname+"*_epc.fits")) + \
                           glob.glob(os.path.join(self.basedir, rootname+"*_spt.fits")) + \
                           glob.glob(os.path.join(self.basedir, rootname+"*_asn.fits")) + \
                           glob.glob(os.path.join(self.basedir, rootname+"*_wav.fits"))
                for item in allfiles:
                    shutil.copy(item, self._sci_dir)
    
        # Run stis_cti
        stis_cti.stis_cti(self._sci_dir, self._dark_dir, self._ref_dir, processes, verbose=True)
        
#-----------------------------------------------------------------------------#

    def copy_products(self):
        """
        Copy the intermediate products from STIS CTI products to self.outdir, 
        since often many calibration iterations are needed and we don't want
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
        for rootname in self.rootname:
            if self.rootname[rootname]["detector"] == "ccd":
                interm_prod = glob.glob(os.path.join(self._sci_dir, rootname+"*flc.fits"))[0]
            else:
                interm_prod = glob.glob(os.path.join(self.basedir, rootname+"*flt.fits"))[0]
            shutil.copy(interm_prod, self.outdir)
            
            print("Copied {0} to {1}".format(interm_prod, self.outdir))
            self.rootname[rootname]["flt"] = glob.glob(os.path.join(self.outdir, rootname+"*"))[0]
        
#-----------------------------------------------------------------------------#

    def shift_spectra(self):
        """
        Shift the FLT and FLC datasets to be aligned in XD using sshift.
        """

        import stistools
        print("Aligning datasets...")

        for rootname in self.rootname:
            myshift = self.rootname[rootname]["shift"] 
            flt = self.rootname[rootname]["flt"]
            if myshift == 0:
                self.rootname[rootname]["shift_flt"] = flt
                continue
            
            # Give the shifted file a new name.
            if "flc.fits" in flt:
                outfile = flt.replace("flc.fits", "shift_flc.fits")
            elif "flt.fits" in flt:
                outfile = flt.replace("flt.fits", "shift_flt.fits")
            self.rootname[rootname]["shift_flt"] = outfile
            
            stistools.sshift.shiftimage(flt, outfile, shift=myshift)
        
#-----------------------------------------------------------------------------#

    def combine_spectra(self, sdqflags=31743):
        """
        Combine the shifted FLTs and FLCs for each visit.
        
        Args:
            sdqflags (int): Bit value descriving serious data quality flags.
        """
        
        print("Combining datasets...")
        import jo_crrej
        from pyraf import iraf
        from iraf import artdata
        from iraf import stsdas, toolbox, imgtools, mstools, mscombine, dqbits
        # Below vv is necessary to use mscombine.
        iraf.reset(use_new_imt="no")
    
        # Turn off DQ=1024, which is not in STIS default SDQFLAGs.
        # This eventually needs to be optimized to use input sdqflags.
        dqbits(bit1=True, bit2=True, bit3=True, bit4=True, bit5=True, bit6=True,
               bit7=True, bit8=True, bit9=True, bit10=True, bit11=False,
               bit12=True, bit13=True, bit14=True, bit15=True, bit16=True)
    
        done_visits = []
        self.combined = {}
        for rootname in self.rootname:
            visit = self.rootname[rootname]["visit"]
            detector = self.rootname[rootname]["detector"]
            # This describes the detector/visit combo. 
            combo = detector + "_" + visit
            
            if combo not in done_visits:
                # Run the custom CR-rejection on CCD data which also 
                # combines images.
                if detector == "ccd":
                    print("#### RUNNING JO'S CR REJECTION!! ####")
                    visit_exposures = [self.rootname[rootname]["shift_flt"] for rootname in self.rootname if self.rootname[rootname]["visit"] == visit if self.rootname[rootname]["detector"] == detector]
                    outfile = os.path.join(self.outdir, "visit{0}_{1}_comb_crc.fits".format(visit, detector))
                    self.rootname[rootname]["comb"] = outfile             
                    print(visit_exposures)
                    jo_crrej.cr_rejection(visit_exposures, outfile, crsigmas=[5], 
                                          initguess="minimum", rej_rad=None, neighbor_thresh=0.8,
                                          sdqflags=sdqflags, update_dq=False, scale_noise=None)  

                # Run mscombine to combine MAMA data.
                else:
                    print("#### RUNNING MSCOMBINE!! ####")
                    visit_exposures = ",".join(self.rootname[rootname]["shift_flt"] for rootname in self.rootname if self.rootname[rootname]["visit"] == visit if self.rootname[rootname]["detector"] == detector)
                    outfile = os.path.join(self.outdir, "visit{0}_{1}_comb_flt.fits".format(visit, detector))
                    mscombine(input=visit_exposures, output=outfile)
                self.combined[outfile] = {"visit": visit, "detector": detector}
    
                done_visits.append(combo)
            
#-----------------------------------------------------------------------------#

    def extract_spectra(self):
        """
        Extract the spectra using x1d.
        """

        from stistools import x1d
        print("Extracting spectra...")
        

        self.x1d = []
        for comb in self.combined: 
            for targ in ["red", "blue"]:
                visit = self.combined[comb]["visit"]
                detector = self.combined[comb]["detector"]
                outfile = os.path.join(self.outdir, "visit{0}_{1}_{2}_x1d.fits".format(
                                       visit, detector, targ))

                # Look up extraction parameters given the observation config.    
                dict_i = self.extract_info[visit][detector][targ]
                x1d.x1d(comb, 
                        output = outfile, 
                        a2center = dict_i["b_spec"], 
                        maxsrch = dict_i["maxsrch"],
                        extrsize = dict_i["height"],
                        bk1offst = dict_i["b_bkg1"] - dict_i["b_spec"],
                        bk2offst = dict_i["b_bkg2"] - dict_i["b_spec"],
                        bk1size = dict_i["b_hgt1"],
                        bk2size = dict_i["b_hgt2"],
                        ctecorr="omit",
                        verbose=True)
            
            print("Extracted spectra in {0}".format(self.outdir))
            self.x1d.append(outfile)
