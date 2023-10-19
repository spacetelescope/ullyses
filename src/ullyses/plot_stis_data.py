import numpy as np
import os
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
try:
    plt.style.use("niceplot.mplstyle")
except OSError:
    pass
from matplotlib import gridspec

COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
          "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"]


def plot_all_x1d(newx1dfile, oldx1dfile, targ, outdir, newname="Custom", oldname="MAST", savepng=False):
    newx1d = fits.getdata(newx1dfile)
    oldx1d = fits.getdata(oldx1dfile)
    grating = fits.getval(newx1dfile, "opt_elem")
    rootname = fits.getval(newx1dfile, "rootname")
    fig0 = compare_dq(newx1d, oldx1d, targ, grating, outdir, newname, oldname, savepng)
    fig1 = overplot(newx1d, oldx1d, targ, grating, outdir, newname, oldname, savepng)
    fig2 = plotdiv(newx1d, oldx1d, targ, grating, outdir, newname, oldname, savepng)
    fig3 = plotdiff(newx1d, oldx1d, targ, grating, outdir, newname, oldname, savepng)
    if savepng is False:
        pdffile = os.path.join(outdir, f"{targ}_{grating}_{rootname.lower()}_1D_compare.pdf")
        pdf = PdfPages(pdffile)
        pdf.savefig(fig0)
        pdf.savefig(fig1)
        pdf.savefig(fig2)
        pdf.savefig(fig3)
        pdf.close()
        print(f"Wrote {pdffile}")
        
        return pdffile
    else:
        return None


def plot_all_2d(flt, acq, x1dfile, targ, outdir, savepng=False):
    spectral_im = fits.getdata(flt)
    rootname = fits.getval(flt, "rootname")
    acq_im = fits.getdata(acq, 4)
    x1d = fits.getdata(x1dfile)
    grating = fits.getval(flt, "opt_elem")
    fig0 = twod_images(spectral_im, acq_im, x1d, targ, grating, outdir, savepng=savepng) 
    fig1 = plot_ee(spectral_im, x1d, targ, grating, outdir, savepng=savepng)
    if savepng is False:
        pdffile = os.path.join(outdir, f"{targ}_{grating}_{rootname.lower()}_2D.pdf")
        pdf = PdfPages(pdffile)
        pdf.savefig(fig0)
        pdf.savefig(fig1)
        pdf.close()
        print(f"Wrote {pdffile}")
        
        return pdffile
    else:
        return None
         

def twod_images(spectral_im, acq_im, x1d, targ, grating, outdir, vmin=-1, 
                vmax=15, log=False, savepng=False):

    if grating == "G230L":
        default = 5
    else:
        default = 3
    extrlocy = x1d["extrlocy"][0]
    center_beg = extrlocy[0]
    center_mid = extrlocy[512]
    center_end = extrlocy[-1]

    fig = plt.figure(figsize=(30, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    if log is True:
        if vmin < 0:
            vmin = 0.1
        vmax = np.log10(vmax)
        spectral_im = np.log10(spectral_im)
        spectral_im = np.nan_to_num(spectral_im, neginf=0)
        acq_im = np.log10(acq_im)
        acq_im = np.nan_to_num(acq_im, neginf=0)
        lab = "log(Counts)"
    else:
        lab = "Counts"
    implot = ax0.imshow(spectral_im, origin="lower", vmin=vmin, vmax=vmax, cmap="inferno", aspect="auto")
    ax0.plot([1, 1023], [center_beg+default, center_end+default], ls="--", color="cyan")
    ax0.plot([1, 1023], [center_beg-default, center_end-default], ls="--", color="cyan")
#    ax0.plot(extrlocy + default, color="cyan")
#    ax0.plot(extrlocy - default, color="cyan")
    ax1.imshow(acq_im, origin="lower", vmin=vmin, vmax=vmax, cmap="inferno")
    fig.colorbar(implot, label=lab)
    
    ax0.set_ylim(center_mid-75, center_mid+75)
    ax0.set_title("Spectral Image With Default Box")
    ax1.set_title("Acq")
    ax0.set_xlabel("X")
    ax0.set_ylabel("Y")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    
    plt.suptitle(f"{targ}")
    
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_2d.png")
        plt.tight_layout()
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plot_ee(spectral_im, x1d, targ, grating, outdir, cols=[250, 500, 750], 
            savepng=False):
    fig,axes0 = plt.subplots(1, 3, figsize=(30,10))
    axes = axes0.flatten()
    for i,col in enumerate(cols):
        data = spectral_im[:, col]
        axes[i].plot(data, color="k")
        axes[i].set_title(f"Column {col}")
        axes[i].set_xlabel("Row")
        axes[i].set_ylabel("Counts")
        
        yloc = int(x1d["extrlocy"][0][col])
        inf = 15
        lo = yloc - inf
        hi = yloc + inf
        total = np.sum(data[lo:hi])
        for j,frac in enumerate([0.85, 0.90, 0.99]):
            ee = total * frac
            running = data[yloc]
            row = 1
            while running < ee:
                running += data[yloc+row]
                running += data[yloc-row]
                row += 1
            axes[i].axvline(yloc+row, color=COLORS[j], linestyle="dashed", 
                            label=f"{frac*100}% EE at height={row*2+1}")
            axes[i].axvline(yloc-row, color=COLORS[j], linestyle="dashed") 

        axes[i].axvline(yloc+inf, color=COLORS[j+1], linestyle="dashed", 
                        label=f"100% EE at height={inf*2+1}")
        axes[i].axvline(yloc-inf, color=COLORS[j+1], linestyle="dashed") 

        if grating == "G230L":
            default = 5
        else:
            default = 3
        axes[i].axvline(yloc+default, color=COLORS[j+2], linestyle="solid", 
                        label=f"Default height={default*2+1}")
        axes[i].axvline(yloc-default, color=COLORS[j+2], linestyle="solid") 
        axes[i].legend(loc="upper right")
        axes[i].set_xlim(yloc-20, yloc+20)
        axes[i].set_ylim(bottom=-20)
    
    plt.suptitle(f"{targ}: Encircled Energies")
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_ee.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def compare_dq(newx1d, oldx1d, targ, grating, outdir, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    newsdq = np.where((newx1d["dq"] & 31743) != 0)
    oldsdq = np.where((oldx1d["dq"] & 31743) != 0)
    ax.plot(oldx1d["wavelength"][0], oldx1d["flux"][0], "k", label=oldname)
    ax.plot(oldx1d["wavelength"][oldsdq], oldx1d["flux"][oldsdq], "rx", markersize=6)
    ax.plot(newx1d["wavelength"][0], newx1d["flux"][0]+1e-14, "royalblue", label=newname)
    ax.plot(newx1d["wavelength"][newsdq], newx1d["flux"][newsdq]+1e-14, "rx", label="SDQ", markersize=6)
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(f"{targ}: Offset Data")
    ax.legend(loc="upper right")
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_dq.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def overplot(newx1d, oldx1d, targ, grating, outdir, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    ax.plot(oldx1d["wavelength"][0], oldx1d["flux"][0], "k", label=oldname)
    ax.plot(newx1d["wavelength"][0], newx1d["flux"][0], "royalblue", label=newname, alpha=0.8)
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(f"{targ}: Overplotted Data")
    ax.legend(loc="upper right")
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_overplot.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plotdiv(newx1d, oldx1d, targ, grating, outdir, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    div = newx1d["flux"][0] / oldx1d["flux"][0]
    ax.plot(newx1d["wavelength"][0], div)
    ax.set_xlabel("Wavelength [A]")
    ax.set_title(f"{targ}: {newname} / {oldname}")
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_div.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plotdiff(newx1d, oldx1d, targ, grating, outdir, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    diff = newx1d["flux"][0] - oldx1d["flux"][0]
    ax.plot(newx1d["wavelength"][0], diff)
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(f"{targ}: {newname} - {oldname}")
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_diff.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig
