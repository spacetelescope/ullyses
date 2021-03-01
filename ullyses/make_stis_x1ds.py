import numpy as np
import os
import glob
from astropy.io import fits as pf
import matplotlib
import matplotlib.pyplot as pl


from calibrate_stis_data import Stisdata

#tts = ["CVSO-104", "CVSO-107", "CVSO-109", "CVSO-146", "CVSO-165", "CVSO-17",
#       "CVSO-176", "CVSO-36", "CVSO-58", "CVSO-90", "V-TX-ORI", "V505-ORI",
#       "V510-ORI"]

tts = ["CVSO-104", "CVSO-107", "CVSO-109", "CVSO-146", "CVSO-17",
       "CVSO-176", "CVSO-36", "CVSO-58", "V-TX-ORI", "V505-ORI",
       "V510-ORI"]

datadir = "/astro/ullyses/tts_dr2"
outdir0 = "out1"

config_dir = "/user/jotaylor/git/ullyses_dp/high_level_science_products/high_level_science_products/config_files"

bad = []
def make_x1ds():
    for targ in tts:
        raws = glob.glob(os.path.join(datadir, targ, "o*_raw.fits"))
        for item in raws:
            hdr = pf.getheader(item)
            grating = hdr["opt_elem"]
            targname = hdr["targname"]
            if targname == targ and grating in ["G750L", "G430L"]:
                infile = item
                outdir = os.path.join(os.path.dirname(item), outdir0)
                config = os.path.join(config_dir, f"{targ.lower()}_{grating.lower()}.yaml")
                if not os.path.exists(config):
                    print(f"NO CONFIG {config} for {targ}")
                    bad.append(item)
                else:
                    print(item, grating, config)
                    S = Stisdata(item, yamlfile=config, outdir=outdir, overwrite=True)
                    S.run_all()
                    this = input("Press enter to continue to next target\n")
    print(f"Files that did not calibrate properly: {bad}")

def check_x1d(newfile, oldfile):
    new = pf.getdata(newfile)
    old = pf.getdata(oldfile)
    compare_dq(new, old)
    overplot(new, old)
    plotdiv(new, old)
    plotdiff(new, old)

def compare_dq(new, old):
    fig,ax = pl.subplots(figsize=(20,7))
    newsdq = np.where((new["dq"] & 31743) != 0)
    oldsdq = np.where((old["dq"] & 31743) != 0)
    ax.plot(old["wavelength"][0], old["flux"][0], "k", label="MAST")
    ax.plot(old["wavelength"][oldsdq], old["flux"][oldsdq], "rx")
    ax.plot(new["wavelength"][0], new["flux"][0]+1e-14, "b", label="New")
    ax.plot(new["wavelength"][newsdq], new["flux"][newsdq]+1e-14, "rx", label="SDQ")
    ax.legend(loc="upper right")
    this = input("")
    pl.close()

def overplot(new, old):
    fig,ax = pl.subplots(figsize=(20,7))
    ax.plot(old["wavelength"][0], old["flux"][0], "k", label="MAST")
    ax.plot(new["wavelength"][0], new["flux"][0], "b", label="New", alpha=0.7)
    ax.legend(loc="upper right")
    this = input("")
    pl.close()

def plotdiv(new, old):
    fig,ax = pl.subplots(figsize=(20,7))
    div = new["flux"][0] / old["flux"][0]
    ax.plot(new["wavelength"][0], div)
    ax.set_title("New - Old")
    this = input("")
    pl.close()

def plotdiff(new, old):
    fig,ax = pl.subplots(figsize=(20,7))
    diff = new["flux"][0] - old["flux"][0]
    ax.plot(new["wavelength"][0], diff)
    ax.set_title("New / Old")
    this = input("")
    pl.close()

if __name__ == "__main__":
    make_x1ds()
