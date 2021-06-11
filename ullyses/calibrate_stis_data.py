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
from stistools import x1d
from stistools.ocrreject import ocrreject
#import stis_cti
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from readwrite_yaml import read_config, write_config
from plot_stis_x1ds import plot_all
os.environ["oref"] = "/grp/hst/cdbs/oref/"

OREF_DIR = "/grp/hst/cdbs/oref"
SYM = "~"
NCOLS = 72
SEP = f"\n!{'~'*70}!\n"


class StisData():
    """
    A class to describe STIS data that need to be custom-calibrated.

    Attributes:
        outdir (str): Path where output products should be written.
        infiles (list): List of input raw datasets.
        nfiles (int): Number of input raw datasets.
        basedir (str): Path to raw STIS data.
        rootname (dict): Dictionary where each key is a dictionary- each key
            is the rootname and vals describe dataset properties (detector,
            dither offset, actual file name).
        _sci_dir (str): 'sci' directory used by stis_cti
        _dark_dir (str): 'dark' directory used by stis_cti
        _ref_dir (str): 'ref' directory used by stis_cti
#!!!        combined (dict): Nested dictionary where keys are the filenames of 
            the combined FLTs and FLCs; each nested dict listing the
            visit and detector values.
        x1d (list): List of final x1d products. 
        
    """
    def __init__(self, infile, yamlfile, outdir=None, overwrite=True): 
        """
        Args:
            infiles: List or wildcard describing all input files.
            outdir: Output directory for products.
        """
        self.raw = self.flt = self.x1d = self.x1d_mast = None
    
        self.do_flag_negatives = self.negatives_flagged = False
        self.custom_dq16_applied = False
        
        self.yamlfile = yamlfile
        config = read_config(yamlfile)
        self.force_dq16 = config["force_dq16"]
        self.target_dict = config["targets"]
        
        self.infile = infile
        self.basedir = os.path.dirname(self.infile)
        if outdir is None:
            nowdt = datetime.datetime.now()
            outdir = os.path.join(self.basedir, nowdt.strftime("%Y%m%d_%H%M"))
        self.outdir = outdir
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        
        self.rootname = fits.getval(infile, "rootname")
        self.visit = self.rootname[4:6]
        self.opt_elem = fits.getval(infile, "opt_elem")
        self.darkfile = fits.getval(infile, "darkfile")
        

    def flag_negatives(self, change_sci_val=False, sci_val=10000, dq=4, thresh=-100):
        """
        For pixels with values below the input negative threshold, change the DQ
        value to a specified value. Also, if specified, change large negative 
        values to large positive values.
    
        Args:
            change_sci_val (Bool): True if large negative pixel SCI values should
                be changed to large positive values, default is False.
            sci_val (Integer or Float): If change_sci_val is True, set large 
                negative values to this value.
            dq (int): DQ value that large negative pixel values should
                be assigned, default is 4.
            thresh (int or float): Threshold below wich pixels should
                be flagged with specified DQ flag, default is -100.
        Returns:
            None
        """

        print("\n", f" FLAGGING NEGATIVES ".center(NCOLS, SYM), "\n")
        print(f"Flagging pixels with counts below {thresh} with DQ={dq}...")
        with fits.open(self.infile, mode="update") as sci_hdu: 
            sci_data = sci_hdu[1].data
            neg = np.where(sci_data <= thresh)
            sci_hdu[3].data[neg] += dq
            if change_sci_val:
                print("Changing large negative values to large positive values...")
                sci_hdu[1].data[neg] = sci_val
        

    def custom_dq16(self, dq=16, threshold=0.06):
        """
        Create custom superdarks that are identical to their CRDS counterparts
        except that they have correct dark-flagging done (STIS reference 
        pipeline does not do it correctly). Apply dark DQ (by default, DQ=16: 
        "Pixel having dark rate >5σ times the median dark level" -STIS IHB) 
        to FLTs.

        Args:
            dq (int): DQ to apply to pixels with large dark values, 
                default is 16.
            threshold (float): If fraction of detector flagged as DQ=16
                is greater than this, perform custom flagging.
        Returns:
            None 
        """

        print("\n", f" CHECKING DQ=16 FLAGS ".center(NCOLS, SYM), "\n")
        
        try:
            shutil.copy(self.infile, self.outdir)
        except shutil.SameFileError:
            pass

        # Read in science FLT dataset.
        with fits.open(self.infile, mode="update") as sci_hdu:
        
            for sciext in range(len(sci_hdu)):
                if sci_hdu[sciext].name != "SCI":
                    continue
                dqext = sciext + 2
                sci_dq = sci_hdu[dqext].data
                darkfile0 = sci_hdu[0].header["darkfile"]
    
                # Check how many pixels are flagged.
                # DQ=512 is "bad pixel in reference file"
                sci_dq16 = np.where((sci_dq&dq == dq) & (sci_dq&512 == 0))
                n_flagged = len(sci_dq16[0])
                total = len(sci_dq) * len(sci_dq[0])
                perc_flagged = n_flagged / total
                if perc_flagged <= 0.06 and self.force_dq16 is False:
                    print(f"For ext={dqext}, less than 6% of pixels are flagged with DQ=16 ({perc_flagged*100.:.2f})")
                    print("Not performing custom DQ flagging")
                    continue
                elif perc_flagged <= 0.06 and self.force_dq16 is True:
                    print(f"For ext={dqext}, less than 6% of pixels are flagged with DQ=16 ({perc_flagged*100.:.2f})")
                    print("But force_dq16 flag is True, so custom DQ flagging will be applied")
                    self.custom_dq16_applied = True
                else:
                    print(f"For ext={dqext}, more than 6% of pixels ({perc_flagged*100.:.2f}) flagged with DQ=16")
                    print("Custom DQ flagging will be applied")
                    self.custom_dq16_applied= True
    
                customdark_dir = os.path.join(self.basedir, "custom_darks")
                if not os.path.isdir(customdark_dir):
                    os.mkdir(customdark_dir)
                    print(f"Made directory: {customdark_dir}")
    
                # Determine DARKFILE filename.
                if self.do_perform_cti is True:
                    darkname = os.path.basename(darkfile0)
                    darkfile = os.path.join(self._ref_dir, darkname)
                else:
                    if "oref" in darkfile0:
                        darkname = darkfile0.split("$")[1]
                        darkfile = os.path.join(OREF_DIR, darkname)
                    else:
                        darkname = os.path.basename(darkfile0)
                        darkfile = os.path.join(OREF_DIR, darkname)


                # Remove any existing DQ={dq} flags from FLT, since they are
                # inherently wrong.
                new_dark_file = os.path.join(customdark_dir, darkname)
                sci_hdu[0].header["DARKFILE"] = new_dark_file
                self.darkfile = new_dark_file
                sci_hdu[dqext].data[sci_dq16] -= dq
                
                written_darks = [os.path.basename(x) for x in 
                                 glob.glob(os.path.join(customdark_dir, "*fits"))]
                # If the custom darkfile has already been written 
                # (from a previous iteration or another dataset close in time), 
                # read from the darkfile and apply the DQ={dq} flags to the FLT. 
                if darkname in written_darks:
                    dark_dq = fits.getdata(new_dark_file, 3)
                    dark_dq16 = np.where(dark_dq&dq == dq)
                    sci_hdu[dqext].data[dark_dq16] |= dq
                    print(f"Already wrote custom superdark:\n{new_dark_file}")
                
                # Create custom darkfile, flag high dark values with dq={dq} 
                else:
                    print(f"Creating superdark with custom DQ={dq} flags")
                    with fits.open(darkfile) as dark_hdu:
                        dark = dark_hdu[1].data
                        dark_dq = dark_hdu[3].data
                        dark_dq16 = np.where(dark_dq&dq == dq)
                        dark_dq[dark_dq16] -= dq
                
                        # Determine 5*stddev + median of the darkfile, anything above this is DQ={dq}. 
                        dark_thresh = (np.std(dark)*5) + np.median(dark)
                        dark_inds = np.where((dark > dark_thresh) | (dark < -dark_thresh))
                        dark_dq[dark_inds] |= dq
                        sci_hdu[dqext].data[dark_inds] |= dq
                
                        new_dark = fits.HDUList( fits.PrimaryHDU() )
                        new_dark[0].header = dark_hdu[0].header
                        for ext in [1,2]:
                            new_dark.append(fits.ImageHDU(dark_hdu[ext].data, dark_hdu[ext].header, name=dark_hdu[ext].name))
                        new_dark.append(fits.ImageHDU(dark_dq, dark_hdu[3].header, name=dark_hdu[3].name))
                        new_dark.writeto(new_dark_file)
                    print(f"Wrote new darkfile: {new_dark_file}")
            

    def extract_spectra(self):
        """
        Extract the spectra using stistools.x1d.
        """

        # Look up extraction parameters for the specific target.
        for target,target_pars in self.target_dict.items():
            print("\n", f" EXTRACTING {target} ".center(NCOLS, SYM), "\n")
            pars = target_pars["x1d"]
        
            outfile = os.path.join(self.outdir, f"{self.rootname}_{target}_x1d.fits")
            self.target_dict[target]["out_x1d"] = outfile
            kwargs = {}
            if "out_drj" not in target_pars:
                if self.detector == "CCD":
                    infile = self.crj
                    kwargs = {"ctecorr": "omit" if self.do_perform_cti is True else "perform"}
                else:
                    infile = self.flt
            else:
                infile = target_pars["out_drj"]
            
            if os.path.exists(outfile):
                os.remove(outfile)
            
            x1d.x1d(infile, 
                output = outfile, 
                a2center = pars["yloc"],
                xoffset = pars["xoffset"], 
                maxsrch = pars["maxsrch"],
                extrsize = pars["height"],
                bk1offst = (pars["b_bkg1"] - pars["yloc"]) if pars["b_bkg1"] is not None else None,
                bk2offst = (pars["b_bkg2"] - pars["yloc"]) if pars["b_bkg2"] is not None else None,
                bk1size = pars["b_hgt1"],
                bk2size = pars["b_hgt2"],
                verbose=True,
                **kwargs)
            print(f"Wrote x1d file: {outfile}")


    def find_product(self, ext):
        prod = os.path.join(self.basedir, self.rootname+"_"+ext+".fits")
        if not os.path.exists(prod):
            prod = os.path.join(self.basedir, "mast_products", self.rootname+"_"+ext+".fits")
            if not os.path.exists(prod):
                prod = None
#            prod = os.path.join(self.outdir, self.rootname+"_"+ext+".fits")
        if prod is not None:
            with fits.open(prod, mode="update") as hdulist:
                hdr0 = hdulist[0].header
                darkfile = hdr0["darkfile"]
                if "ctiref" in darkfile:
                    darkfile = darkfile.replace("$ctirefs/", self._ref_dir)
                    hdr0.set("DARKFILE", darkfile)
        return prod


    def update_header(self):
        nowdt = datetime.datetime.now()
        nowdt_str = nowdt.strftime("%Y%m%d_%H%M")
        for target,target_pars in self.target_dict.items():
            with fits.open(target_pars["out_x1d"], mode="update") as hdulist:
                hdulist[0].header["HLSP_LVL"] = 0
                hdulist[0].header["ULL_CAL"] = "COMPLETE"
                hdulist[0].header["CAL_DATE"] = nowdt_str
                hdulist[0].header["TARGNAME"] = target.upper()
                hdulist[0].header["FILENAME"] = os.path.basename(target_pars["out_x1d"])
                if "coords" in target_pars:
                    hdulist[0].header["RA_TARG"] = target_pars["coords"]["ra"]
                    hdulist[0].header["DEC_TARG"] = target_pars["coords"]["dec"]


    def printintro(self):
        print("\n", f" RUNNING ON {os.path.basename(self.infile)} ".center(NCOLS, SYM), "\n")
        print(f"Filename: {self.infile}")
        print(f"Detector: {self.detector}")
        print(f"Grating: {self.opt_elem}")
        print("Target(s):")
        for target in self.target_dict:
            print(f"\t{target}")


    def compare_x1ds(self):
        if self.x1d_mast is None:
            return
        print("\n", f" CREATING DIAGNOSTIC PLOTS ".center(NCOLS, SYM), "\n")
        mastx1d = self.x1d_mast
        for target,target_pars in self.target_dict.items():
            customx1d = target_pars['out_x1d']
            pdffile = plot_all(customx1d, mastx1d, target, self.outdir)
            self.target_dict[target]["out_plot"] = pdffile


class StisCcd(StisData):
    def __init__(self, infile, yamlfile, outdir=None, overwrite=True):
        super().__init__(infile, yamlfile, outdir, overwrite)
        self.detector = "CCD"
        
        prod = os.path.join(self.basedir, self.rootname+"_sx1.fits")
        if os.path.exists(prod):
            self.x1d_mast = prod
        prod = os.path.join(self.basedir, self.rootname+"_flt.fits")
        if os.path.exists(prod):
            self.flt = prod
        self.crj = None
        self.infile_type = "crj"
        if self.infile_type == "crj":
            self.crj = infile
    
        self.do_custom_dq16 = True
        self.do_perform_cti = self.custom_cti_applied = False
        self.cti_proc = 10
        self.needs_crrej = self.do_crrej = self.custom_crrej_applied = False
        self.defringe_applied = False
        
        self._sci_dir = self.outdir
        self._dark_dir = "/astro/ullyses/stis_ccd_data/darks"
        self._ref_dir = "/astro/ullyses/stis_ccd_data/cti_refs"
        os.environ["ctirefs"] = "/astro/ullyses/stis_ccd_data/cti_refs/"
        
        self.crsplit = fits.getval(infile, "crsplit")
   
    
    def perform_cti(self):
        """
        Run the STIS CTI code on STIS CCD data.
        """
        
        print("\n", f" PERFORMING CTI ".center(NCOLS, SYM), "\n")
        
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
            try:
                shutil.copy(item, self._sci_dir)
            except shutil.SameFileError:
                pass
    
        # Run stis_cti
        stis_cti.stis_cti(self._sci_dir, self._dark_dir, self._ref_dir, self.cti_proc, verbose=True, clean=True)
        self.flt = os.path.join(self.outdir, self.rootname+"_flc.fits")
        self.crj = os.path.join(self.outdir, self.rootname+"_crc.fits")
#        self.copy_products()


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
        crc_files = glob.glob(os.path.join(self._sci_dir, "*crc.fits"))
        for item in flc_files+crc_files:
            shutil.copy(item, self.outdir)
        
        print(f"Copied FLC and CRC files to {self.outdir}")
        self.flt = os.path.join(self.outdir, self.rootname+"_flc.fits")
        self.crj = os.path.join(self.outdir, self.rootname+"_crc.fits")


    def check_crrej(self):
        
        print("\n", f" CHECKING CR REJECTION ".center(NCOLS, SYM), "\n")
        
        if self.x1d_mast is None:
            print("WARNING: No x1d files found, cannot check CR rejection rate")
            return
        if self.flt is None:
            print("WARNING: No FLT files found, cannot check CR rejection rate")
            return

        x1d_data = fits.getdata(self.x1d_mast)[0]
        extr_mask = np.zeros((1024, 1024))
        del_pix = x1d_data["EXTRSIZE"]/2.
        for column in range(1024):
            row_mid =  x1d_data['EXTRLOCY'][column] - 1 #EXTRLOCY is 1-indexed.
            gd_row_low = int(np.ceil( row_mid - del_pix))
            gd_row_high = int(np.floor(row_mid + del_pix))
            extr_mask[gd_row_low:gd_row_high+1,column] = 1
        n_tot = np.count_nonzero(extr_mask) * self.crsplit

        n_rej = []
        with fits.open(self.flt) as flt_hdu:
            for i in range(3,len(flt_hdu)+1 ,3):
                flt_data = flt_hdu[i].data
                rej = flt_data[ (extr_mask == 1) & (flt_data & 8192 != 0)] #Data quality flag 8192 (2^13) used for CR rejected pixels
                n_rej.append(np.count_nonzero(rej))

        t_exp = float(fits.getval(self.x1d_mast, "texptime"))
        rej_rate = float(fits.getval(self.x1d_mast, "rej_rate"))

        # Calculate the rejection fraction and the rate of rejected pixels per sec
        n_pix_rej = np.sum(np.array(n_rej))
        frac_rej = n_pix_rej/float(n_tot)
        rej_rate = n_pix_rej * 1024.0**2 / (t_exp * n_tot)

        # Calculate the expected rejection fraction rate given the rate in header
        # This is the rej_rate * ratio of number of pixels in full CCD vs extract region (n_tot /CRSPLIT)
        frac_rej_expec =  rej_rate * t_exp /(1024.*1024.*self.crsplit)

        print('Percentage of pixels rejected as CRs')          
        print(f'   Extraction Region: {frac_rej*100:.2f}%')
        print(f'   Full CCD Frame: {frac_rej_expec*100:.2f}% \n')

        if frac_rej < (frac_rej_expec * 5) or frac_rej < .05:
            print("Extraction region rejection rate is not high enough to trigger custom CR rejection")
            return
        else:
            print("Extraction region rejection rate is above threshold")
            self.needs_crrej = True

    
    def custom_crrej(self):
        """
        Check on CR rejection rate comes from Joleen Carlberg's code: 
        https://github.com/spacetelescope/stistools/blob/jkc_cr_analysis/stistools/crrej_exam.py
        """
       
        print("\n", f" PERFORMING CR REJECTION ".center(NCOLS, SYM), "\n")

        if self.do_crrej is True and self.needs_crrej is True and self.infile_type == "flt":
            print("WARNING!!! CUSTOM CRREJ NOT YET IMPLEMENTED")
            pass
        else:
            return
        
        if self.do_perform_cti is True:
            outfile = os.path.join(self.outdir, self.rootname+"_crc.fits")
        else:
            outfile = os.path.join(self.outdir, self.rootname+"_crj.fits")
        if os.path.exists(outfile):
            os.remove(outfile)
        
        # Look up crrej parameters given the target config file.
        ocrreject(self.flt,
                  output=outfile,
                  initgues=self.crrej_pars["initgues"],
                  crsigmas=self.crrej_pars["crsigmas"],
                  crradius=self.crrej_pars["crradius"],
                  crthresh=self.crrej_pars["crthresh"],
                  crmask=self.crrej_pars["crmask"],
                  verbose=self.crrej_pars["verbose"])

        print(f"Wrote crj file: {outfile}")
        self.crj = outfile

    
    def defringe(self):
        
        for target,target_pars in self.target_dict.items():
            print("\n", f" DEFRINGING {target} ".center(NCOLS, SYM), "\n")
            if target_pars["do_defringe"] is False:
                return
            if self.opt_elem not in ["G750L", "G750M"]:
                return
            self.do_defringe = True 
            
            pars = target_pars["defringe"]
            fringeflat = pars["fringeflat"]
            rawfringe = os.path.join(self.basedir, fringeflat)
            fringeroot = os.path.basename(fringeflat)[:9]
            
            outnorm = os.path.join(self.outdir, fringeroot+"_nsp.fits")
            with fits.open(rawfringe, mode="update") as hdulist:
                hdr0 = hdulist[0].header
                darkfile = fits.getval(self.infile, "darkfile")
                hdr0.set("DARKFILE", darkfile)
            if os.path.exists(outnorm):
                os.remove(outnorm)
            
            stistools.defringe.normspflat(inflat=rawfringe,
                       do_cal=True,
                       outflat=outnorm)
    
            if self.opt_elem == "G750L":
                with fits.open(outnorm, mode="update") as hdulist:
                    hdulist[1].data[:,:250] = 1
    
            outmk = os.path.join(self.outdir, fringeroot+"_mff.fits")
            if os.path.exists(outmk):
                os.remove(outmk)
            stistools.defringe.mkfringeflat(inspec=self.crj,
                         inflat=outnorm,
                         outflat=outmk,
                         do_shift=pars["do_shift"],
                         beg_shift=pars["beg_shift"],
                         end_shift=pars["end_shift"],
                         shift_step=pars["shift_step"],
                         do_scale=pars["do_scale"],
                         beg_scale=pars["beg_scale"],
                         end_scale=pars["end_scale"],
                         scale_step=pars["scale_step"],
                         opti_spreg=pars["opti_spreg"],
                         rms_region=pars["rms_region"],
                         extrloc=pars["extrloc"],
                         extrsize=pars["extrsize"])

            outfile = stistools.defringe.defringe(science_file=self.crj,
                               fringe_flat=outmk,
                               overwrite=True,
                               verbose=True)
            outfile_name = outfile.replace("_drj", f"_{target}_drj")
            outfile_dest = os.path.join(self.outdir, os.path.basename(outfile_name))
            self.target_dict[target]["out_drj"] = outfile_dest

            try:
                shutil.copyfile(outfile, outfile_dest)
            except shutil.SameFileError:
                pass
            print(f"Wrote defringed crj file: {outfile_dest}")


    def run_all(self):
        self.printintro()
        if self.do_perform_cti is True:
            self.perform_cti()
        if self.do_flag_negatives is True:
            self.flag_negatives()
        self.check_crrej()
        if self.do_crrej is True:
            self.custom_crrej()
        if self.do_custom_dq16 is True:
            self.custom_dq16()
        self.defringe()
        self.extract_spectra()
        self.update_header()
        self.compare_x1ds()
        self.printfinal()
    

    def printfinal(self):
        print("\n", f" RECALIBRATION SUMMARY ".center(NCOLS, SYM), "\n")
        print("Target(s):")
        for target in self.target_dict:
            print(f"\t{target}")
        print(f"YAML file: {self.yamlfile}")
        print(f"Input file: {self.infile}")
        print(f"FLT: {self.flt}")
        print(f"CRJ: {self.crj}")
        print(f"DRJ(s):")
        for target,target_pars in self.target_dict.items():
            if "out_drj" in target_pars:
                print(f"\t{target_pars['out_drj']}")
        print(f"Custom X1D file(s):")
        for target,target_pars in self.target_dict.items():
            print(f"\t{target_pars['out_x1d']}")
        print(f"MAST X1D file: {self.x1d_mast}")
        print(f"Comparison plots:")
        for target,target_pars in self.target_dict.items():
            print(f"\t{target_pars['out_plot']}")

        print("")
        if self.do_perform_cti is True:
            print("* Pixel-based CTI correction performed")
        else:
            print("* Default CalSTIS empirical CTI correction performed")
        if self.custom_dq16_applied is True:
            print("* Custom DQ=16 flagging has been applied")
        else:
            print("* Default DQ=16 flagging was used")
        if self.needs_crrej is True:
            print("* WARNING!!! CR rejection rate is higher than allowed")
        else:
            print("* Default CR rejection parameters were used")
        if self.do_defringe is False:
            print("* No fringe correction applied")
        else:
            print("* Fringe correction applied")
        print("* Custom extraction performed")


class StisMama(StisData):
    def __init__(self, infile, yamlfile, outdir=None, overwrite=True):
        super().__init__(infile, yamlfile, outdir, overwrite)
        self.detector = "MAMA"
        self.do_custom_dq16 = True
        self.do_defringe = False
        prod = os.path.join(self.basedir, self.rootname+"_x1d.fits")
        if os.path.exists(prod):
            self.x1d_mast = prod
        self.infile_type = "flt"
        if self.infile_type == "flt":
            self.flt = infile


    def run_all(self):
        self.printintro()
        if self.do_flag_negatives is True:
            self.flag_negatives()
        if self.do_custom_dq16 is True:
            self.custom_dq16()
        self.extract_spectra()
        self.update_header()
        self.compare_x1ds()
        self.printfinal()
    

    def printfinal(self):
        print("\n", f" RECALIBRATION SUMMARY ".center(NCOLS, SYM), "\n")
        print("Target(s):")
        for target in self.target_dict:
            print(f"\t{target}")
        print(f"YAML file: {self.yamlfile}")
        print(f"Input file: {self.infile}")
        print(f"FLT(s): {self.flt}")
        print(f"Custom X1D file(s):")
        for target,target_pars in self.target_dict.items():
            print(f"\t{target_pars['out_x1d']}")
        print(f"MAST X1D file: {self.x1d_mast}")
        print(f"Comparison plots:")
        for target,target_pars in self.target_dict.items():
            print(f"\t{target_pars['out_plot']}")

        print("")
        if self.custom_dq16_applied is True:
            print("* Custom DQ=16 flagging has been applied")
        else:
            print("* Default DQ=16 flagging was used")
        print("* Custom extraction performed")


def wrapper(indir, yamlfile, outdir=None, overwrite=True):
    config = read_config(yamlfile)
    infile = os.path.join(indir, config["infile"])
    detector = fits.getval(infile, "detector")
    opt_elem = fits.getval(infile, "opt_elem")
    if detector == "CCD" and opt_elem !="MIRVIS":
        S = StisCcd(infile, yamlfile, outdir, overwrite)
    else:
        S = StisMama(infile, yamlfile, outdir, overwrite)
    S.run_all()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", 
                        help="Path to input STIS dataset")
    parser.add_argument("-y", "--yaml", 
                        help="Name of YAML configuration file")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Directory for output products")
    parser.add_argument("-c", "--clobber", default=False,
                        action="store_true",
                         help="If True, overwrite existing products")
    args = parser.parse_args()
    wrapper(args.indir, args.yaml, args.outdir, args.clobber)

