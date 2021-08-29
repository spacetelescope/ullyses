from astropy.io import fits
import calcos
import os
import shutil
import glob
import pandas as pd
import datetime

SHIFTS = {"sz10": "cos_shift_files/sz10_shifts.txt"}
CUSTOM_DIR = "/astro/ullyses/custom_cal"
VETTED_DIR = "/astro/ullyses/all_vetted_data_dr3"
OUTDIR_ROOT = None
nowdt = datetime.datetime.now()
if OUTDIR_ROOT is None:
    OUTDIR_ROOT = nowdt.strftime("%Y%m%d_%H%M")

def apply_shifts():
    for targ,shift_file in SHIFTS.items():
        seen = []
        df = pd.read_csv(shift_file, delim_whitespace=True, 
            names=["rootname", "fppos", "flash", "segment", "shift1"], 
            header=None)
        roots = df["rootname"]
        for root in roots:
            ipppss = root[:6]
            if ipppss not in seen:
                seen.append(ipppss)
            else:
                continue
            files = glob.glob(os.path.join(CUSTOM_DIR, targ,  
                          f"{ipppss}*asn.fits"))
            assert len(files) != 0, f"No files found for {targ}"
            outdir = os.path.join(CUSTOM_DIR, targ, OUTDIR_ROOT)
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            calcos.calcos(files[0], shift_file=shift_file, outdir=outdir, 
                          verbosity=0)
            x1ds = glob.glob(os.path.join(outdir, "*x1d.fits"))
            copydir = os.path.join(VETTED_DIR, targ)
            for x1d in x1ds:
                existing = os.path.join(copydir, os.path.basename(x1d))
                if os.path.exists(existing):
                    os.remove(existing)
                with fits.open(x1d, mode="update") as hdulist:
                    hdulist[0].header["HLSP_LVL"] = 0
                shutil.copy(x1d, copydir)
            print(f"Copied custom x1ds to {copydir}")
    print(f"Copied custom x1ds to {VETTED_DIR}")

if __name__ == "__main__":
    apply_shifts()
