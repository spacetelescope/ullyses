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

COLORS = ["#9970ab", "#5aae61", "#d95f02", "#e7298a", "#a6cee3", "#1f78b4", 
          "#fb9a99", "#fdbf6f", "#ffe44b", 
          "#b15928", "#cab2d6", "#b2df8a", "#000000", "#7a7a7a", "#911cff", 
          "#1bddf2"]
CONTRAST_COLORS = ["#5a0a7d" ,"#0d4f13"]
SDQFLAGS = 31743 


def plot_all_x1d(newx1dfile, oldx1dfile, targ, outdir, newname="Custom", oldname="MAST", savepng=False):
    newx1d = fits.getdata(newx1dfile)
    oldx1d = fits.getdata(oldx1dfile)
    grating = fits.getval(newx1dfile, "opt_elem")
    rootname = fits.getval(newx1dfile, "rootname")
    filename = os.path.basename(oldx1dfile)
    fig0 = compare_dq(newx1d, oldx1d, targ, grating, outdir, filename, newname, oldname, savepng)
    fig1 = overplot(newx1d, oldx1d, targ, grating, outdir, filename, newname, oldname, savepng)
    fig2 = plotdiv(newx1d, oldx1d, targ, grating, outdir, filename, newname, oldname, savepng)
    fig3 = plotdiff(newx1d, oldx1d, targ, grating, outdir, filename, newname, oldname, savepng)
    if savepng is False:
        pdffile = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_{rootname.lower()}_1d_compare.pdf")
        pdf = PdfPages(pdffile)
        pdf.savefig(fig0)
        pdf.savefig(fig1)
        pdf.savefig(fig2)
        pdf.savefig(fig3)
        pdf.close()
        print(f"Wrote {pdffile}")
        fig0.clear()
        fig1.clear()
        fig2.clear()
        fig3.clear()
        plt.close()         
        return pdffile
    else:
        fig0.clear()
        fig1.clear()
        fig2.clear()
        fig3.clear()
        plt.close()         
        return None


def plot_all_2d(flt, acq, x1dfile, targ, outdir, savepng=False):
    rootname = fits.getval(flt, "rootname")
    x1d = fits.getdata(x1dfile)
    grating = fits.getval(flt, "opt_elem")
    fig0 = twod_images(flt, acq, x1d, targ, grating, outdir, savepng=savepng) 
    fig1 = plot_ee(flt, x1d, targ, grating, outdir, savepng=savepng)
    if savepng is False:
        pdffile = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_{rootname.lower()}_2d.pdf")
        pdf = PdfPages(pdffile)
        pdf.savefig(fig0)
        pdf.savefig(fig1)
        pdf.close()
        print(f"Wrote {pdffile}")
        
        return pdffile
    else:
        return None
         

def twod_images(flt, acq, x1d, targ, grating, outdir, 
                log=False, savepng=False):
    spectral_im = fits.getdata(flt)
    if acq is not None:
        acq_im = fits.getdata(acq, 4)
        vmax_acq = np.median(acq_im) + np.median(acq_im)*1
    vmin = -1
    vmax_im = 15
    detector = fits.getval(flt, "detector") 
    if detector != "CCD":
        default = 5.5
    else:
        default = 3.5
    extrlocy = x1d["extrlocy"][0]
    center_beg = extrlocy[0] - 1 # make it 0 indexed 
    center_mid = extrlocy[512] - 1
    center_end = extrlocy[-1] - 1

    fig = plt.figure(figsize=(30, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    if log is True:
        if vmin < 0:
            vmin = 0.1
        vmax_im = np.log10(vmax_im)
        if acq is not None:
            vmax_acq = np.log10(vmax_acq)
            acq_im = np.log10(acq_im)
            acq_im = np.nan_to_num(acq_im, neginf=0)
        spectral_im = np.log10(spectral_im)
        spectral_im = np.nan_to_num(spectral_im, neginf=0)
        lab = "log(Counts)"
    else:
        lab = "Counts"
    implot = ax0.imshow(spectral_im, origin="lower", vmin=vmin, vmax=vmax_im, cmap="inferno", aspect="auto")
    ax0.plot(extrlocy + default - 1, color="cyan", linestyle="dashed")
    ax0.plot(extrlocy - default - 1, color="cyan", linestyle="dashed")
    if acq is not None:
        ax1.imshow(acq_im, origin="lower", vmin=vmin, vmax=vmax_acq, cmap="inferno")
        ax1.set_title(f"Acquisition Image\n{os.path.basename(acq)}")
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
#    fig.colorbar(implot, label=lab)
    
    ax0.set_ylim(center_mid-75, center_mid+75)
    ax0.set_title(f"Spectral Image With Default Box: {os.path.basename(flt)}")
    ax0.set_xlabel("X")
    ax0.set_ylabel("Y")
    
    plt.suptitle(f"{targ} {grating}")
    plt.tight_layout()
    
    if savepng is True:
        png = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_2d.png")
        plt.tight_layout()
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plot_ee(flt, x1d, targ, grating, outdir, cols=[250, 500, 750], 
            savepng=False):
    fig,axes0 = plt.subplots(1, 3, figsize=(30,10))
    axes = axes0.flatten()
    spectral_im = fits.getdata(flt)
    for i,col in enumerate(cols):
        data = spectral_im[:, col]
        axes[i].plot(data, color="k")
        axes[i].set_title(f"X = {col}")
        axes[i].set_xlabel("Y (Row)")
        axes[i].set_ylabel("Counts")
        
        yloc = int(x1d["extrlocy"][0][col])
        inf = 15
        lo = yloc - inf
        hi = yloc + inf
        total = np.sum(data[lo:hi])
        for j,frac in enumerate([0.85, 0.90, 0.99]):
            ee = total * frac
            try:
                running = data[yloc]
            except IndexError:
                continue
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

        detector = fits.getval(flt, "detector") 
        if detector != "CCD":
            default = 5.5
        else:
            default = 3.5
        axes[i].axvline(yloc+default, color="grey", linestyle="solid", 
                        label=f"Default height={default*2}")
        axes[i].axvline(yloc-default, color="grey", linestyle="solid") 
        axes[i].legend(loc="upper right")
        axes[i].set_xlim(yloc-20, yloc+20)
        axes[i].set_ylim(bottom=-20)
    
    plt.suptitle(f"{targ} {grating}: Encircled Energies")
    plt.tight_layout()
    if savepng is True:
        png = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_ee.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plot_one_x1d(x1dfile, targ, outdir, savefig=False):
    x1d = fits.getdata(x1dfile)
    grating = fits.getval(x1dfile, "opt_elem")
    rootname = fits.getval(x1dfile, "rootname")
    filename = os.path.basename(x1dfile)
    grating = fits.getval(x1dfile, "opt_elem")
    fig,ax = plt.subplots(figsize=(20,7))
    sdq = np.where((x1d["dq"] & SDQFLAGS) != 0)
    ax.plot(x1d["wavelength"][0], x1d["flux"][0], COLORS[0])
    ax.plot(x1d["wavelength"][sdq], x1d["flux"][sdq], "x", 
            color=CONTRAST_COLORS[0], label="SDQ", markersize=6)
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(f"{targ} {grating}: {filename}")
    ax.legend(loc="upper right")
    if savefig is True:
        png = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_{rootname.lower()}_1d.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        fig.clear()
        plt.close()
        return png
    fig.clear()
    plt.close()         

def compare_dq(newx1d, oldx1d, targ, grating, outdir, filename, newname="Custom", 
               oldname="MAST", savepng=False):
    fig,axes = plt.subplots(2, 1, figsize=(20,14), sharex=True)
    ax0, ax1 = axes.flatten()
    newsdq = np.where((newx1d["dq"] & SDQFLAGS) != 0)
    oldsdq = np.where((oldx1d["dq"] & SDQFLAGS) != 0)
    offset = np.median(oldx1d["flux"][0])*3
    for ax in [ax0, ax1]: 
        ax.plot(oldx1d["wavelength"][0], oldx1d["flux"][0], COLORS[1], 
                label=oldname)
        ax.plot(oldx1d["wavelength"][oldsdq], oldx1d["flux"][oldsdq], "x", 
                color=CONTRAST_COLORS[1], markersize=6)
        ax.plot(newx1d["wavelength"][0], newx1d["flux"][0]+offset, COLORS[0], 
                label=newname)
        ax.plot(newx1d["wavelength"][newsdq], newx1d["flux"][newsdq]+offset, "x", 
                color=CONTRAST_COLORS[0], label="SDQ", markersize=6)
    ylo = np.median(oldx1d["flux"][0]) - np.median(oldx1d["flux"][0])*2
    yhi = np.median(newx1d["flux"][0]) + np.median(newx1d["flux"][0])*2 + offset
    ax1.set_ylim(ylo, yhi)
    ax1.set_xlabel("Wavelength [A]")
    ax0.set_ylabel("Flux")
    ax0.set_title(f"{targ} {grating}: {filename} Offset Data")
    ax1.set_title("Zoom in")
    ax0.legend(loc="upper right")
    ax1.legend(loc="upper right")
    plt.tight_layout()
    if savepng is True:
        png = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_dq.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def overplot(newx1d, oldx1d, targ, grating, outdir, filename, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    ax.plot(oldx1d["wavelength"][0], oldx1d["flux"][0], COLORS[1], 
            label=oldname, alpha=.8)
    ax.plot(newx1d["wavelength"][0], newx1d["flux"][0], COLORS[0], 
            label=newname, alpha=0.8)
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(f"{targ} {grating}: {filename} Overplotted Data")
    ax.legend(loc="upper right")
    plt.tight_layout()
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_overplot.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plotdiv(newx1d, oldx1d, targ, grating, outdir, filename, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    div = newx1d["flux"][0] / oldx1d["flux"][0]
    ax.plot(newx1d["wavelength"][0], div, COLORS[0])
    ax.set_xlabel("Wavelength [A]")
    ax.set_title(f"{targ} {grating} {filename}: {newname} / {oldname}")
    plt.tight_layout()
    if savepng is True:
        png = os.path.join(outdir, f"{targ}_{grating}_div.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig 


def plotdiff(newx1d, oldx1d, targ, grating, outdir, filename, newname="Custom", oldname="MAST", savepng=False):
    fig,ax = plt.subplots(figsize=(20,7))
    diff = newx1d["flux"][0] - oldx1d["flux"][0]
    ax.plot(newx1d["wavelength"][0], diff, COLORS[0])
    ax.set_xlabel("Wavelength [A]")
    ax.set_ylabel("Flux")
    ax.set_title(f"{targ} {grating} {filename}: {newname} - {oldname}")
    plt.tight_layout()
    if savepng is True:
        png = os.path.join(outdir, f"{targ.lower()}_{grating.lower()}_diff.png")
        plt.savefig(png, bbox_inches="tight")
        print(f"Wrote {png}")
        plt.cla()
        plt.close()
    return fig


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--x1d", type=str,
                        help="Name of 1D spectrum (x1d or sx1)")
    parser.add_argument("--flt", type=str,
                        help="Name of 2D spectral image (flt or crj)")
    parser.add_argument("--acq", type=str,
                        help="Name of acquisition image")
    parser.add_argument("--targ", type=str,
                        help="Name of target")
    parser.add_argument("--outdir", type=str,
                        help="Name of directory to save plots to")
    args = parser.parse_args()
    plot_one_x1d(args.x1d, args.targ, args.outdir, savefig=True)
    plot_all_2d(args.flt, args.acq, args.x1d, args.targ, args.outdir, savepng=False)
