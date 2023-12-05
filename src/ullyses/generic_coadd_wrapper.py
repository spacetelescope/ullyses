from astropy.io import fits
from collections import defaultdict
import os

from ullyses.coadd import COSSegmentList, STISSegmentList, FUSESegmentList, CCDSegmentList  

def coadd_files(infiles, outdir, outfile=None, clobber=False):
     
    nonvofiles = [x for x in infiles if "_vo.fits" not in x]
    vofiles = [x for x in infiles if "_vo.fits" in x]
    # collect the gratings that we will loop through
    # coadd.py will find the correct files itself,
    # but we need to know which gratings are present
    uniqmodes = defaultdict(list)
    
    for infile in nonvofiles:
        prihdr = fits.getheader(infile)
        obsmode = (prihdr['INSTRUME'], prihdr['OPT_ELEM'], prihdr['DETECTOR'])
        uniqmodes[obsmode].append(infile)

    if vofiles:
        if len(vofiles) != 1:
            print("More than 1 FUSE data file, aborting")
        else:
            obsmode = ('FUSE', 'FUSE', 'FUSE')
            uniqmodes[obsmode].append(vofiles[0])

    if not uniqmodes:
        print(f'No data to coadd, EXITING')
        return

    level = 2
    for obskey in uniqmodes:
        instrument, grating, detector = obskey
        infiles = uniqmodes[obskey]
        # this instantiates the class
        if instrument == 'COS':
            prod = COSSegmentList(grating, infiles=infiles)
        elif instrument == 'STIS':
            if detector == 'CCD':
                prod = CCDSegmentList(grating, infiles=infiles)
            else:
                prod = STISSegmentList(grating, infiles=infiles)
        elif instrument == 'FUSE':
            prod = FUSESegmentList(grating, infiles=infiles)
        else:
            print(f'Unknown mode [{instrument}, {grating}, {detector}]')
            return
        if len(prod.members) > 0:
            if outdir is None:
                outdir = indir
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            if outfile is None:
                grating_outfile = f"{instrument.lower()}_{grating.lower()}_coadd.fits"
                grating_outfile = os.path.join(outdir, grating_outfile)
            prod.create_output_wavelength_grid()
            prod.coadd()
            prod.target = prod.get_targname()
            prod.targ_ra, prod.targ_dec, prod.coord_epoch = prod.get_coords()
            prod.write(grating_outfile, clobber)
            print(f"Wrote {grating_outfile}")
