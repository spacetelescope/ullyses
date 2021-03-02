import numpy as np
import os
import glob
from astropy.io import fits as pf
import matplotlib
import matplotlib.pyplot as pl
from stistools import x1d

from calibrate_stis_data import Stisdata

#tts = ["CVSO-104", "CVSO-107", "CVSO-109", "CVSO-146", "CVSO-165", "CVSO-17",
#       "CVSO-176", "CVSO-36", "CVSO-58", "CVSO-90", "V-TX-ORI", "V505-ORI",
#       "V510-ORI"]

tts = ["CVSO-104", "CVSO-107", "CVSO-109", "CVSO-146", "CVSO-17",
       "CVSO-176", "CVSO-36", "CVSO-58", "CVSO-90", "V-TX-ORI", "V505-ORI",
       "V510-ORI"]

datadir = "/astro/ullyses/tts_dr2"
outdir0 = "20210301"

config_dir = "/user/jotaylor/git/ullyses_dp/high_level_science_products/high_level_science_products/config_files"

bad = []
def make_ccd_x1ds():
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
                    print("")
                    print("#"*80)
                    if S.sx1_mast is not None:
                        check_x1d(S.sx1, S.sx1_mast, f"{targ}_{grating}", outdir)
#                    else:
#                        check_x1d(S.sx1, S.sx1, targ, outdir)
                    print("#"*80)
                    print("")
    print(f"Files that did not calibrate properly: {bad}")

def make_mama_x1ds():
    outfile = os.path.join(datadir, "CVSO-104", outdir0, "oe9k1s020_nonsci410_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-104/mast_products/oe9k1s020_flt.fits"),
            output=outfile,
            a2center=410,
            maxsrch=10,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "CVSO-109", outdir0, "oe9k2s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-109/mast_products/oe9k2s010_flt.fits"),
            output=outfile,
            a2center=502,
            maxsrch=0)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    outfile = os.path.join(datadir, "CVSO-109", outdir0, "oe9k2s010_nonsci516_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-109/mast_products/oe9k2s010_flt.fits"),
            output=outfile,
            a2center=516.3,
            maxsrch=0,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "CVSO-36", outdir0, "oe9k5s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-36/mast_products/oe9k5s010_flt.fits"),
            output=outfile,
            a2center=502.2578,
            maxsrch=0,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "CVSO-17", outdir0, "oe9k3s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-36/mast_products/oe9k3s010_flt.fits"),
            output=outfile,
            a2center=502.73044,
            maxsrch=0,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "CVSO-58", outdir0, "oe9j3s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-58/mast_products/oe9j3s010_flt.fits"),
            output=outfile,
            extrsize=9)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "CVSO-165", outdir0, "oe9j2s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-165/mast_products/oe9j2s010_flt.fits"),
            output=outfile,
            a2center=514.5,
            maxsrch=0,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    outfile = os.path.join(datadir, "CVSO-165", outdir0, "oe9j2s010_nonsci527_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-165/mast_products/oe9j2s010_flt.fits"),
            output=outfile,
            a2center=527.0,
            maxsrch=0,
            extrsize=9)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    outfile = os.path.join(datadir, "CVSO-165", outdir0, "oe9j2s010_nonsci733_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-165/mast_products/oe9j2s010_flt.fits"),
            output=outfile,
            a2center=733,
            maxsrch=13,
            extrsize=7,
            bk1offst=-521,
            bk2offst=79)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "CVSO-176", outdir0, "oe9k4s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-176/mast_products/oe9k4s010_flt.fits"),
            output=outfile,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    
    outfile = os.path.join(datadir, "V505-ORI", outdir0, "oe9i3s010_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "V505-ORI/mast_products/oe9i3s010_flt.fits"),
            output=outfile,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")

def rename_targs():
    comps = {"CVSO-14": {"GAIA-DR3-3217634157789741952": (83.02660499856, -1.18336721131)}, 
        "CVSO-109": {"CVSO-109B": (83.13599552466, -1.22957460015)}, 
#        "CVSO-165": {"GAIA-DR3-3217473697810165504": (084.75978669326, -1.34115742897),
#                     "CVSO-165B": (,)},
        "CVSO-36": {"CVSO-36B": (81.45900692210, 1.82725967070)}}
     
    for targ in comps:
        d = comps[targ]
        for comp in d:
            files = glob.glob(os.path.join(datadir, targ, outdir0, f"*{comp}*x1d.fits"))
            for item in files:
                with pf.open(item, mode="update") as hdulist:
                    hdulist[0].header["RA_TARG"] = d[comp][0]
                    hdulist[0].header["DEC_TARG"] = d[comp][1]

    mains = {"CVSO-165": [("oe9j2s010_x1d.fits", "CVSO-165A"),
                          ("oe9j2s020_x1d.fits", "CVSO-165A"),
                          ("oe9j2s030_x1d.fits", "CVSO-165A")],
             "CVSO-36": [("oe9k5s010_x1d.fits", "CVSO-36A"),
                         ("oe9k5s020_x1d.fits", "CVSO-36A"),
                         ("oe9k5s030_x1d.fits", "CVSO-36A")],
             "CVSO-109": [("oe9k2s010_x1d.fits", "CVSO-109A"),
                          ("oe9k2s020_x1d.fits", "CVSO-109A"),
                          ("oe9k2s030_x1d.fits", "CVSO-109A")]}
    for targ in mains:
        d = mains[targ]
        files = glob.glob(os.path.join(datadir, targ, outdir0, d[0]))
        if len(files) > 0:
            print(f"something went wrong with {targ}")
            continue
        x1d = files[0]
        with pf.open(x1d, mode="update") as hdulist:
            hdulist[0].header["TARGNAME"] = d[1]

def check_x1d(newfile, oldfile, targ, outdir):
    new = pf.getdata(newfile)
    old = pf.getdata(oldfile)
    spl = newfile.split("/")
    newname = spl[-2]+"/"+spl[-1]
    spl = oldfile.split("/")
    oldname = spl[-2]+"/"+spl[-1]
    compare_dq(new, newname, old, oldname, targ, outdir)
    overplot(new, newname, old, oldname, targ, outdir)
    plotdiv(new, newname, old, oldname, targ, outdir)
    plotdiff(new, newname, old, oldname, targ, outdir)

def compare_dq(new, newname, old, oldname, targ, outdir):
    fig,ax = pl.subplots(figsize=(20,7))
    newsdq = np.where((new["dq"] & 31743) != 0)
    oldsdq = np.where((old["dq"] & 31743) != 0)
    ax.plot(old["wavelength"][0], old["flux"][0], "k", label=oldname)
    ax.plot(old["wavelength"][oldsdq], old["flux"][oldsdq], "rx")
    ax.plot(new["wavelength"][0], new["flux"][0]+1e-14, "royalblue", label=newname)
    ax.plot(new["wavelength"][newsdq], new["flux"][newsdq]+1e-14, "rx", label="SDQ")
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(targ)
    ax.legend(loc="upper right")
    png = os.path.join(outdir, f"{targ}_dq.png")
    pl.savefig(png, bbox_inches="tight")
    print(f"Wrote {png}")
    pl.cla()
    pl.close()

def overplot(new, newname, old, oldname, targ, outdir):
    fig,ax = pl.subplots(figsize=(20,7))
    ax.plot(old["wavelength"][0], old["flux"][0], "k", label=oldname)
    ax.plot(new["wavelength"][0], new["flux"][0], "royalblue", label=newname, alpha=0.8)
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(targ)
    ax.legend(loc="upper right")
    png = os.path.join(outdir, f"{targ}_overplot.png")
    pl.savefig(png, bbox_inches="tight")
    print(f"Wrote {png}")
    pl.cla()
    pl.close()

def plotdiv(new, newname, old, oldname, targ, outdir):
    fig,ax = pl.subplots(figsize=(20,7))
    div = new["flux"][0] / old["flux"][0]
    ax.plot(new["wavelength"][0], div)
    ax.set_xlabel("Wavelength [A]")
    ax.set_title(f"{targ}: New - Old")
    png = os.path.join(outdir, f"{targ}_div.png")
    pl.savefig(png, bbox_inches="tight")
    print(f"Wrote {png}")
    pl.cla()
    pl.close()

def plotdiff(new, newname, old, oldname, targ, outdir):
    fig,ax = pl.subplots(figsize=(20,7))
    diff = new["flux"][0] - old["flux"][0]
    ax.plot(new["wavelength"][0], diff)
    ax.set_xlabel("Wavelength [A]")
    ax.set_title(f"{targ}: New / Old")
    png = os.path.join(outdir, f"{targ}_diff.png")
    pl.savefig(png, bbox_inches="tight")
    print(f"Wrote {png}")
    pl.cla()
    pl.close()

if __name__ == "__main__":
    make_ccd_x1ds()
    make_mama_x1ds()
    rename_targs()
