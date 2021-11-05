import numpy as np
from astropy.io import fits
import glob
import os
import shutil
import calcos

import splittag_wrapper
import timeseries

CODEDIR = os.path.dirname(__file__)
VERSION = "dr4"
ALLDATA_DIR = "/astro/ullyses/ULLYSES_DATA"
VETTED_DIR = "/astro/ullyses/all_vetted_data_dr4"
CUSTOM_DIR = f"/astro/ullyses/custom_cal/{VERSION}"

TARGS = ["v-tw-hya", "v-bp-tau", "v-ru-lup", "v-gm-aur"]

BINS = {"v-tw-hya": {"g160m": {"time": 30, "wave": 3, "min_exptime": 20}, # exptime = 300
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}, # exptime = 30
        "v-ru-lup": {"g160m": {"time": 50, "wave": 3, "min_exptime": 40}, # exptime = 220
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}, # exptime = 30
        "v-bp-tau": {"g160m": {"time": 60, "wave": 3, "min_exptime": 40}, # exptime = 128
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}, # exptime = 196
        "v-gm-aur": {"g160m": {"time": 90, "wave": 6, "min_exptime": 50}, # exptime = 186
                     "g230l": {"time": 10, "wave": 1, "min_exptime": 9}}} # exptime = 184

BAD_IPPPSS = ["le9d1k", # TW Hydra
              "leit1d", "leitad", "leit1l", # RU Lup
              "lek71f"] # GM Aur 


G230L_DISPTAB = os.path.join(CODEDIR, 'ctts_recalibration/ullyses_cos_nuv_disp.fits')

WL_SHIFT = {'le9d1c': "ctts_recalibration/twhya_shifts.txt", 
            'le9d1g': "ctts_recalibration/twhya_shifts.txt"} 
# TW Hydra exposures
#['le9d1cdeq',
# 'le9d1cdeq',
# 'le9d1cdgq',
# 'le9d1cdgq',
# 'le9d1cdiq',
# 'le9d1cdiq',
# 'le9d1cdkq',
# 'le9d1cdkq',
# 'le9d1cdmq',
# 'le9d1cdmq',
# 'le9d1cdoq',
# 'le9d1cdoq',
# 'le9d1gw7q']

def copy_origdata(targs=TARGS):
    """
    Copy the original raw data, to ensure nothing is mistakenly edited.
    Data is copied from /astro/ullyses/ULLYSES_DATA/{targ} to
    /astro/ullyses/custom_cal/{version}/{targ}.

    Args:
        targs (list): List of targets for which to copy data.
    """

    for targ in targs:
        targ_u = targ.upper()
        datadir = os.path.join(ALLDATA_DIR, targ_u)
        files = glob.glob(os.path.join(datadir, "l*corrtag*fits"))
        files += glob.glob(os.path.join(datadir, "l*rawtag*.fits"))
        files += glob.glob(os.path.join(datadir, "l*spt*.fits"))
        files += glob.glob(os.path.join(datadir, "l*asn*.fits"))
        files += glob.glob(os.path.join(datadir, "l*x1d.fits"))
        destdir = os.path.join(CUSTOM_DIR, targ)
        if not os.path.exists(destdir):
            os.makedirs(destdir)
        for item in files:
            shutil.copy(item, destdir)
        print(f"Copied original files to {destdir}")

def calibrate_data(targs=TARGS):
    """
    Calibrate data which require special calibration. This will be for
    1) all G230L data which need to be calibrated with Will's custom DISPTAB
    2) some TW Hydra exposures which were offset in wavelength. 
    """

    calrequired0 = []
    for targ in targs:
        datadir = os.path.join(CUSTOM_DIR, targ)
        raws = glob.glob(os.path.join(datadir, "l*rawtag*fits"))
        for raw in raws:
            if fits.getval(raw, "opt_elem") == "G230L":
                with fits.open(raw, mode="update") as hdulist:
                    hdr0 = hdulist[0].header
                    hdr0.set("DISPTAB", G230L_DISPTAB)
                root = os.path.basename(raw)[:9]
                calrequired0.append(root)
        wl_shift_files = list(WL_SHIFT.keys())
        calrequired0 += wl_shift_files
        calrequired = [x[:6] for x in calrequired0]

        asns = glob.glob(os.path.join(datadir, "*asn.fits"))
        outdir = os.path.join(CUSTOM_DIR, targ, "calcos_wl")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        for asn in asns:
            ipppss = os.path.basename(asn)[:6]
            if ipppss not in calrequired:
                continue
            if ipppss in BAD_IPPPSS: # This is a bad visit
                continue
            if ipppss in WL_SHIFT: 
                calcos.calcos(asn, shift_file=WL_SHIFT[ipppss], outdir=outdir, verbosity=0)
            else:
                calcos.calcos(asn, outdir=outdir, verbosity=0)

        print(f"\n Calibrated data that required a shift or custom DISPTAB\n, output dir={outdir}")


def copy_caldata(targs=TARGS):
    """
    Copy the products for each target.
    """

    for targ in targs:
        datadir = os.path.join(CUSTOM_DIR, targ)
        # First copy the calibrated output
        outdir = os.path.join(CUSTOM_DIR, targ, "calcos_wl")
        corrs = glob.glob(os.path.join(outdir, "*corrtag*"))
        copydir = os.path.join(CUSTOM_DIR, targ, VERSION)
        copydirfuv = os.path.join(copydir, "g160m")
        copydirnuv = os.path.join(copydir, "g230l")
        copydirfuvx1d = os.path.join(copydirfuv, "exp")
        copydirnuvx1d = os.path.join(copydirnuv, "exp")
        for d in [copydirfuv, copydirnuv, copydirfuvx1d, copydirnuvx1d]:
            if not os.path.exists(d):
                os.makedirs(d)
        for item in corrs:
            if fits.getval(item, "opt_elem") == "G160M":
                d = copydirfuv
            else:
                d = copydirnuv
#            filename = os.path.basename(item)
#            sptfile = filename.split("_")[0]+"_spt.fits"
#            spt = os.path.join(datadir, sptfile)
            shutil.copy(item, d)
#            shutil.copy(spt, d)
            
            root = os.path.basename(item)[:9]
            x1ds = glob.glob(os.path.join(outdir, root+"*x1d.fits"))
            for x1d in x1ds:
                shutil.copy(x1d, os.path.join(d, "exp"))

        # Then copy the corrtags that didn't require custom calibration
        orig_corrs = glob.glob(os.path.join(datadir, "*corrtag*"))
        orig_corrfiles = [os.path.basename(x) for x in orig_corrs]
        corrfiles = [os.path.basename(x) for x in corrs]
        for i in range(len(orig_corrs)):
            ipppss = os.path.basename(orig_corrs[i])[:6]
            if orig_corrfiles[i] not in corrfiles and ipppss not in BAD_IPPPSS:
                if fits.getval(orig_corrs[i], "opt_elem") == "G160M":
                    d = copydirfuv
                else:
                    d = copydirnuv
                shutil.copy(orig_corrs[i], d)

        x1ds = glob.glob(os.path.join(outdir, "*x1d.fits"))
        orig_x1ds = glob.glob(os.path.join(datadir, "*x1d.fits"))
        orig_x1dfiles = [os.path.basename(x) for x in orig_x1ds]
        x1dfiles = [os.path.basename(x) for x in x1ds]
        for i in range(len(orig_x1ds)):
            ipppss = os.path.basename(orig_x1ds[i])[:6]
            if orig_x1dfiles[i] not in x1dfiles and ipppss not in BAD_IPPPSS: 
                if fits.getval(orig_corrs[i], "opt_elem") == "G160M":
                    d = copydirfuv
                else:
                    d = copydirnuv
                shutil.copy(orig_corrs[i], os.path.join(d, "exp"))

        print(f"\nCopied all corrtags and x1ds to {copydir}/g160m/ and g230l/\n")


def create_splittags(targs=TARGS):
    for targ in targs:
        datadir = os.path.join(CUSTOM_DIR, targ)

        for grat in ["g160m", "g230l"]:
            indir = os.path.join(CUSTOM_DIR, targ, VERSION, grat)
            outdir = os.path.join(indir, "split")
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            # First copy SPT files
            infiles = glob.glob(os.path.join(indir, "*corrtag*fits"))
            for infile in infiles:
                root = os.path.basename(infile)[:9]
                spt = root+"_spt.fits"
                copied = os.path.join(outdir, spt)
                if not os.path.exists(copied):
                    sptfile = os.path.join(datadir, spt)
                    shutil.copy(sptfile, outdir)

            t = BINS[targ][grat]["time"]
            splittag_wrapper.main(indir, outdir, incr=t, numcores=10)
            calcos_out = glob.glob(os.path.join(outdir, "calcosout", "*x1d.fits"))
            for item in calcos_out:
                shutil.move(item, outdir)


def correct_vignetting(targs=TARGS):
    for targ in targs:
        indirs = [os.path.join(CUSTOM_DIR, targ, VERSION, "g230l", "exp"),
                  os.path.join(CUSTOM_DIR, targ, VERSION, "g230l", "split")]
        for indir in indirs:
            files = glob.glob(os.path.join(indir, "*x1d.fits"))
            for item in files:
                if fits.getval(item, "cenwave") == 2950:
                    root = fits.getval(item, "rootname").lower()
                    scale_file = os.path.join(CODEDIR, "ctts_recalibration", f"{root}_scale.txt")
                    assert os.path.exists(scale_file), f"No scaling file found for {item}"
                    scale = np.loadtxt(scale_file)
                    with fits.open(item, mode="update") as hdulist:
                        assert len(scale) == len(hdulist[1].data["flux"][1]), f"Shape of FITS and scaling factor do not match for {item}"
                        hdulist[1].data["flux"][1] /= scale # NUVB is 0th index

        print(f'Applied scaling factor to G230L/2950 NUVB data in {os.path.join(CUSTOM_DIR, targ, VERSION, "g230l")}') 


def create_timeseries(targs=TARGS):
    for targ in targs:
        for grat in ["g160m", "g230l"]:
            outdir = f"/astro/ullyses/ULLYSES_HLSP/{targ}/{VERSION.lower()}"
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            indir = os.path.join(CUSTOM_DIR, targ, VERSION, grat, "exp")
            outfile = os.path.join(outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{VERSION.lower()}_tss.fits")
            timeseries.process_files(grat.upper(), outfile, indir, overwrite=True) 
            
            indir = os.path.join(CUSTOM_DIR, targ, VERSION, grat, "split")
            outfile = os.path.join(outdir, f"hlsp_ullyses_hst_cos_{targ}_{grat}_{VERSION.lower()}_split-tss.fits")
            timeseries.process_files(grat.upper(), outfile, indir, 
                    BINS[targ][grat]["wave"], BINS[targ][grat]["min_exptime"], overwrite=True) 


def main():
    copy_origdata()
    calibrate_data()
    copy_caldata()
    create_splittags()
    correct_vignetting()
    create_timeseries()

if __name__ == "__main__":
    main()
