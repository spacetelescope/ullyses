from astropy.io import fits
import numpy as np

"""
This code adds a data quality array to FUSE data, which
does not come with one as standard. Users can specify
wavelength ranges to be marked as bad the add_dq_col function.
"""


def add_column(infile, outfile, colext, colname, colformat, colvals, 
               colunit=None, overwrite=False):
    """
    Add a new column to a FITS file with HDU Tables as extenstions.

    Args:
        infile (str): File to modify
        outfile (str): Name of output modified file
        colext (int): Extension to add column to
        colname (str): Name of new column
        colformat (str): Format of new column
        colvals (array): Array of new column values
        colunit (str): Optional, unit of new column 
        overwrite (Bool): If True, overwrite any existing files
    """

    hdr0 = fits.getheader(infile, 0)
    new_hdu0 = fits.PrimaryHDU(header=hdr0)

    with fits.open(infile) as tmp:
        nextend = len(tmp)
    new_tables = []
    for ext in range(1, nextend):
        data = fits.getdata(infile, ext)
        cols = data.columns
        names = data.names
        formats = data.formats
        units = [x.unit for x in cols]
        new_cols = []
        for i in range(len(cols)):
            c = fits.Column(name=names[i], format=formats[i], 
                    array=data[names[i]], unit=units[i])
            new_cols.append(c)
        if ext == colext:
            c = fits.Column(name=colname, format=colformat,
                    array=colvals, unit=colunit)
            new_cols.append(c)
        coldefs = fits.ColDefs(new_cols)

        hdr = fits.getheader(infile, ext)
        t = fits.BinTableHDU.from_columns(coldefs, header=hdr)
        new_tables.append(t)

    final_hdus = [new_hdu0] + new_tables
    new_hdulist = fits.HDUList(final_hdus)
    new_hdulist.writeto(outfile, overwrite=overwrite)
    print(f"Wrote {outfile}")


def add_dq_col(infile, outfile, wlstart, wlend, dqflag, overwrite=False):
    """
    Add a DQ column for FUSE data.

    Args: 
        infile (str): File to modify
        outfile (str): Name of output modified file
        wlstart (int, float, or array-like): Beginning of wavelength range(s)
            to flag as bad DQ
        wlend (int, float, or array-like): End of wavelength range(s)
            to flag as bad DQ
        dqflag (int, float, or array-like): DQ flag(s) corresponding to each
            wavelength range to flag.
        overwrite (Bool): If True, overwrite any existing files
    """
    if not isinstance(wlstart, (list, np.ndarray)):
        wlstart = [wlstart]
    if not isinstance(wlend, (list, np.ndarray)):
        wlend = [wlend]
    if not isinstance(dqflag, (list, np.ndarray)):
        dqflag = [dqflag]
    wlstart += [0, 1179.9]
    wlend += [912, -1]
    dqflag += [2, 2]

    good_dq = 0
    
    data = fits.getdata(infile, 1)
    wl = data["wave"][0]
    arrlen = len(wl)
    dqarr = np.zeros(arrlen).astype(int)
    dqarr.fill(good_dq)
    for i in range(len(wlstart)):
        start = wlstart[i]
        end = wlend[i]
        bad_dq = dqflag[i]
        if end == -1:
            end = 99999
        if start == 0:
            start = -99999
        inds = np.where((wl >= start) & (wl <= end))
        dqarr[inds] = bad_dq
    dqarr = dqarr.reshape(1, arrlen)
    add_column(infile, outfile, colext=1, colname="DQ", colformat=str(arrlen)+"I", 
               colvals=dqarr, colunit=None, overwrite=overwrite)
