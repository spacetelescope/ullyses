import shutil
import numpy as np
import os
import glob
from astropy.io import fits as pf
import matplotlib
import matplotlib.pyplot as pl
from stistools import x1d
import subprocess

from coadd import STISSegmentList
from calibrate_stis_data import Stisdata
from fuse_add_dq import add_column

tts = ["CVSO-104", "CVSO-107", "CVSO-109", "CVSO-146", "CVSO-165", "CVSO-17",
       "CVSO-176", "CVSO-36", "CVSO-58", "CVSO-90", "V-TX-ORI", "V505-ORI",
       "V510-ORI"]

version = "dr2"
hlspdir = "/astro/ullyses/ULLYSES_HLSP"
drdir = "/astro/ullyses/all_vetted_data_dr2"
datadir = "/astro/ullyses/tts_dr2"
outdir0 = "v2"

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
    subprocess.run(["chmod", "-R", "777", datadir])
    print(f"Files that did not calibrate properly: {bad}")

def make_mama_x1ds():
    outfile = os.path.join(datadir, "CVSO-104", outdir0, "oe9k1s020_GAIA-DR3-3217634157789741952_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-104/mast_products/oe9k1s020_flt.fits"),
            output=outfile,
            xoffset=2.1295,
            a2center=410,
            maxsrch=10,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    outfile = os.path.join(datadir, "CVSO-104", outdir0, "oe9k1s020_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-104/mast_products/oe9k1s020_flt.fits"),
            output=outfile,
            xoffset=2.1295)
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
    outfile = os.path.join(datadir, "CVSO-109", outdir0, "oe9k2s010_CVSO-109B_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-109/mast_products/oe9k2s010_flt.fits"),
            output=outfile,
            xoffset=22.0115,
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
    x1d.x1d(os.path.join(datadir, "CVSO-17/mast_products/oe9k3s010_flt.fits"),
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
            maxsrch=5,
            extrsize=7)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    outfile = os.path.join(datadir, "CVSO-165", outdir0, "oe9j2s010_CVSO-165B_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-165/mast_products/oe9j2s010_flt.fits"),
            output=outfile,
            a2center=527.0,
            xoffset=-.9214,
            maxsrch=5,
            extrsize=9)
    with pf.open(outfile, mode="update") as hdulist:
        hdulist[0].header["HLSP_LVL"] = 0
    print("#"*80, "\n")
    outfile = os.path.join(datadir, "CVSO-165", outdir0, "oe9j2s010_GAIA-DR3-3217473697810165504_x1d.fits")
    if os.path.exists(outfile):
        os.remove(outfile)
    x1d.x1d(os.path.join(datadir, "CVSO-165/mast_products/oe9j2s010_flt.fits"),
            output=outfile,
            a2center=733,
            maxsrch=7,
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
    subprocess.run(["chmod", "-R", "777", datadir])

def rename_targs():
    comps = {"CVSO-104": {"GAIA-DR3-3217634157789741952": (83.02660499856, -1.18336721131)}, 
        "CVSO-109": {"CVSO-109B": (83.13599552466, -1.22957460015)}, 
        "CVSO-165": {"GAIA-DR3-3217473697810165504": (84.75978669326, -1.34115742897),
                     "CVSO-165B": (84.760675, -1.34227778)},
        "CVSO-36": {"CVSO-36B": (81.45900692210, 1.82725967070)}}
     
    for targ in comps:
        d = comps[targ]
        for comp in d:
            files = glob.glob(os.path.join(datadir, targ, outdir0, f"*{comp}*x1d.fits"))
            for item in files:
                with pf.open(item, mode="update") as hdulist:
                    hdulist[0].header["RA_TARG"] = d[comp][0]
                    hdulist[0].header["DEC_TARG"] = d[comp][1]
                    hdulist[0].header["TARGNAME"] = comp

    mains = {"CVSO-109": [("oe9k2s010_x1d.fits", "CVSO-109A"),
                          ("oe9k2s020_x1d.fits", "CVSO-109A"),
                          ("oe9k2s030_x1d.fits", "CVSO-109A")],
            "CVSO-165": [("oe9j2s010_x1d.fits", "CVSO-165A"),
                          ("oe9j2s020_x1d.fits", "CVSO-165A"),
                          ("oe9j2s030_x1d.fits", "CVSO-165A")]}
    for targ in mains:
        mapping = mains[targ]
        for item in mapping:
            files = glob.glob(os.path.join(datadir, targ, outdir0, item[0]))
            if len(files) > 1:
                print(f"something went wrong with {targ}")
                continue
            if len(files) == 0: #already changed filename
                continue
            x1d = files[0]
            with pf.open(x1d, mode="update") as hdulist:
                hdulist[0].header["TARGNAME"] = item[1]
            newname = x1d.replace("x1d.fits", f"{item[1]}_x1d.fits")
            shutil.move(x1d, newname)
            print(f"Moved {x1d} to {newname}")

    filenames = {"CVSO-165": [("oe9j2s010_CVSO-165_x1d.fits", "oe9j2s010_x1d.fits"),
                      ("oe9j2s020_CVSO-165_x1d.fits", "oe9j2s020_x1d.fits"),
                      ("oe9j2s030_CVSO-165_x1d.fits", "oe9j2s030_x1d.fits"),
                      ("oe9j2s010_CVSO-165A_x1d.fits", "oe9j2s010_x1d.fits"),
                      ("oe9j2s020_CVSO-165A_x1d.fits", "oe9j2s020_x1d.fits"),
                      ("oe9j2s030_CVSO-165A_x1d.fits", "oe9j2s030_x1d.fits"),
                      ("oe9j2s010_CVSO-165B_x1d.fits", "oe9j2s010_x1d.fits"),
                      ("oe9j2s020_CVSO-165B_x1d.fits", "oe9j2s020_x1d.fits"),
                      ("oe9j2s030_CVSO-165B_x1d.fits", "oe9j2s030_x1d.fits"),
                      ("oe9j2s010_GAIA-DR3-3217473697810165504_x1d.fits", "oe9j2s010_x1d.fits"),
                      ("oe9j2s020_GAIA-DR3-3217473697810165504_x1d.fits", "oe9j2s020_x1d.fits"),
                      ("oe9j2s030_GAIA-DR3-3217473697810165504_x1d.fits", "oe9j2s030_x1d.fits")],
             "CVSO-109": [("oe9k2s010_CVSO-109_x1d.fits", "oe9k2s010_x1d.fits"),
                      ("oe9k2s020_CVSO-109_x1d.fits", "oe9k2s020_x1d.fits"),
                      ("oe9k2s030_CVSO-109_x1d.fits", "oe9k2s030_x1d.fits"),
                      ("oe9k2s010_CVSO-109A_x1d.fits", "oe9k2s010_x1d.fits"),
                      ("oe9k2s020_CVSO-109A_x1d.fits", "oe9k2s020_x1d.fits"),
                      ("oe9k2s030_CVSO-109A_x1d.fits", "oe9k2s030_x1d.fits"),
                      ("oe9k2s010_CVSO-109B_x1d.fits", "oe9k2s010_x1d.fits"),
                      ("oe9k2s020_CVSO-109B_x1d.fits", "oe9k2s020_x1d.fits"),
                      ("oe9k2s030_CVSO-109B_x1d.fits", "oe9k2s030_x1d.fits")],
             "CVSO-36": [("oe9k5s020_CVSO-36B_x1d.fits", "oe9k5s020_x1d.fits"),
                      ("oe9k5s030_CVSO-36B_x1d.fits", "oe9k5s030_x1d.fits")],
             "CVSO-104": [("oe9k1s020_GAIA-DR3-3217634157789741952_x1d.fits", "oe9k1s020_x1d.fits"),
                      ("oe9k1s030_GAIA-DR3-3217634157789741952_x1d.fits", "oe9k1s030_x1d.fits"),
                      ("oe9k1s040_GAIA-DR3-3217634157789741952_x1d.fits", "oe9k1s040_x1d.fits")]}

    for targ in filenames:
        fmapping = filenames[targ]
        for item in fmapping:
            files = glob.glob(os.path.join(datadir, targ, outdir0, item[0]))
            if len(files) > 1:
                print(f"something went wrong with {targ}")
                continue
            x1d = files[0]
            key = item[1]
            with pf.open(x1d, mode="update") as hdulist:
                hdulist[0].header["FILENAME"] = key
    subprocess.run(["chmod", "-R", "777", datadir])



def copy_mama_x1ds():
    for targ in tts:
        outdir_files = glob.glob(os.path.join(datadir, targ, outdir0, "*x1d.fits"))
        gratings = [pf.getval(x, "detector") for x in outdir_files]
        if "NUV-MAMA" not in gratings:
            mama_files = glob.glob(os.path.join(datadir, targ, "mast_products", "*010_x1d.fits"))
            if len(mama_files) != 1:
                print(f"something went wrong with MAMA checks for {targ}")
            mama_file = mama_files[0]
            shutil.copy(mama_file, os.path.join(datadir, targ, outdir0))
            print(f"Copied {mama_file} to {os.path.join(datadir, targ, outdir0)}")

def copy_files():
    files = glob.glob(os.path.join(datadir, "*", outdir0, "*x1d.fits"))
    for item in files:
        targ0 = pf.getval(item, "targname")
        targ = targ0.lower()
        dest = os.path.join(drdir, targ)
        if not os.path.exists(dest):
            os.makedirs(dest)
        x1d = os.path.basename(item)
        sx1 = x1d.replace("x1d.fits", "sx1.fits")
        sx1file = os.path.join(dest, sx1)
        if os.path.exists(sx1file):
            os.remove(sx1file)
            print(f"Removed {sx1file}")
        destx1d = os.path.join(dest, x1d)
        if os.path.exists(destx1d):
            os.remove(destx1d)
            print(f"Removed {destx1d}")
        shutil.copy(item, dest)
        print(f"Copied {item} to {dest}")

class STIScoadd(STISSegmentList):
    def create_output_wavelength_grid(self):
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
        self.output_dq = np.zeros(self.nelements).astype(int)
        self.output_gross = np.zeros(self.nelements)
        self.output_net = np.zeros(self.nelements)
        self.output_sumgross = np.zeros(self.nelements)
        self.output_sumnet = np.zeros(self.nelements)
        self.output_sumerrors = np.zeros(self.nelements)
        for i,segment in enumerate(self.members):
            if self.datasets[i] == ignore_dq_file:
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
            if self.datasets[i] != ignore_dq_file:
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
  

def coadd_1d_spectra():
    targs = {"CVSO-109": {"G430L": ("oe9k2s020_CVSO-109A_x1d.fits", "oe9k2s020_CVSO-109B_x1d.fits"),
                          "G750L": ("oe9k2s030_CVSO-109A_x1d.fits", "oe9k2s030_CVSO-109B_x1d.fits"),
                          "G230L": ("oe9k2s010_CVSO-109A_x1d.fits", "oe9k2s010_CVSO-109B_x1d.fits")},
             "CVSO-165": {"G230L": ("oe9j2s010_CVSO-165A_x1d.fits", "oe9j2s010_CVSO-165B_x1d.fits")}}

    for targ in targs:
        d = targs[targ]
        for grating in d:
            files = d[grating]
            filenames = [os.path.join(datadir, targ, outdir0, f) for f in files]
            coadd_dir = os.path.join(datadir, targ, outdir0, f"{grating}_coadd")
            if not os.path.exists(coadd_dir):
                os.mkdir(coadd_dir)
            for item in filenames:
                shutil.copy(item, coadd_dir)
                print(f"Copied {item} to {coadd_dir}")
            root = files[0][:9]
            combined0 = f"{root}_{targ}_x1d.fits"
            combined = os.path.join(datadir, targ, outdir0, combined0)
            if os.path.exists(combined):
                os.remove(combined)
                print(f"Removed {combined}")
            prod = STIScoadd(grating, path=coadd_dir, weighting_method='unity')
            prod.target = prod.ull_targname()
            prod.targ_ra, prod.targ_dec = prod.ull_coords()
            prod.create_output_wavelength_grid()
            ignore_file = os.path.join(coadd_dir, files[1])
            prod.coadd(ignore_dq_file=ignore_file)

            nelements = len(prod.output_net)
            wl_arr = prod.output_wavelength.reshape(1, nelements)
            flux_arr = prod.output_flux.reshape(1, nelements)
            err_arr = prod.output_errors.reshape(1, nelements)
            net_arr = prod.output_net.reshape(1, nelements)
            gross_arr = prod.output_gross.reshape(1, nelements)
            dq_arr = prod.output_dq.reshape(1, nelements)
            hdr0 = pf.getheader(filenames[0], 0)
            hdr0["TARGNAME"] = targ
            new_hdu0 = pf.PrimaryHDU(header=hdr0)
            cols = []
            cols.append(pf.Column(name="WAVELENGTH", format=f"{nelements}D", 
                array=wl_arr, unit="Angstroms"))
            cols.append(pf.Column(name="FLUX", format=f"{nelements}E", 
                array=flux_arr, unit="erg/s/cm**2/Angstrom"))
            cols.append(pf.Column(name="ERROR", format=f"{nelements}E", 
                array=err_arr, unit="erg/s/cm**2/Angstrom"))
            cols.append(pf.Column(name="NET", format=f"{nelements}E", 
                array=net_arr, unit="Counts/s"))
            cols.append(pf.Column(name="GROSS", format=f"{nelements}E", 
                array=gross_arr, unit="Counts/s"))
            cols.append(pf.Column(name="DQ", format=f"{nelements}I", 
                array=dq_arr, unit=None))
            coldefs = pf.ColDefs(cols)
            hdr = pf.getheader(filenames[0], 1)
            t = pf.BinTableHDU.from_columns(coldefs, header=hdr)
            final_hdus = [new_hdu0, t]
            new_hdulist = pf.HDUList(final_hdus)
            new_hdulist.writeto(combined, overwrite=True)
            print(f"Wrote {combined}")
    subprocess.run(["chmod", "-R", "777", datadir])
#            add_column(combined, combined, 1, "WAVELENGTH", f"{nelements}D", wl_arr, colunit="Angstroms", overwrite=True)
#            add_column(combined, combined, 1, "FLUX", f"{nelements}E", flux_arr, colunit="erg/s/cm**2/Angstrom", overwrite=True)
#            add_column(combined, combined, 1, "ERROR", f"{nelements}E", err_arr, colunit="erg/s/cm**2/Angstrom", overwrite=True)
#            add_column(combined, combined, 1, "NET", f"{nelements}E", net_arr, colunit="Counts/s", overwrite=True)
#            add_column(combined, combined, 1, "GROSS", f"{nelements}E", gross_arr, colunit="Counts/s", overwrite=True)
#            add_column(combined, combined, 1, "DQ", f"{nelements}I", dq_arr, overwrite=True)

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

def copy_yamlfiles():
    yamlfiles0 = glob.glob("config_files/*yaml")
    yamlfiles = [x for x in yamlfiles0 if os.path.basename(x) != "target_grating.yaml"]
    for item in yamlfiles:
        f = os.path.basename(item)
        spl = f.split("_")
        targ = spl[0]
        grating = spl[1].split(".")[0]
        newname = f"hlsp_ullyses_hst_stis_{targ}_{grating}_{version}_spec.yaml"
        dest = os.path.join(hlspdir, targ, "dr2")
        shutil.copyfile(item, os.path.join(dest, newname))
        print(f"Copied {item} to {os.path.join(dest, newname)}")


if __name__ == "__main__":
    make_ccd_x1ds()
    make_mama_x1ds()
    copy_mama_x1ds()
    rename_targs()
    coadd_1d_spectra()
    copy_files()
    copy_yamlfiles()
