import argparse
import os
import glob
import numpy as np

from astropy.io import fits

from coadd import COSSegmentList

version = 'v0.1'

'''
Currently coadd.py only works for COS and only combines obs of the same grating.
This wrapper goes through each target folder in the ullyses data directory and find
the COS data and which gratings are present. This info is then feeded into coadd.py.
'''


def main(indir, outdir):

    for root, dirs, files in os.walk(indir, topdown=False):

        print(root)
        targetname = root.split('/')[-1]
        print(f"   {targetname}")

        # collect the gratings that we will loop through
        # coadd.py will find the correct files itself,
        # but we need to know which gratings are present
        modes = []

        for myfile in glob.glob(os.path.join(root, 'l*_x1d.fits')):  # only grabbing COS with the l*_x1d.fits
            f1 = fits.open(myfile)
            prihdr = f1[0].header
            modes.append(prihdr['OPT_ELEM'])

        if not modes:
            print(f'No COS data to coadd for {targetname}.')
            continue

        uniqmodes = np.unique(modes)

        for grating in uniqmodes:

            # this initiates the class
            prod = COSSegmentList(grating, path=root)

            # these two calls perform the main functions
            if len(prod.members) > 0:
                prod.create_output_wavelength_grid()
                prod.coadd()
                prod.target = targetname.lower()
                # this writes the output file
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                outname = create_output_file_name(prod)
                outname = outdir + '/' + outname
                prod.write(outname)
                print(f"   Wrote {outname}")
            else:
                print(f"No valid data for grating {grating}")

def create_output_file_name(prod):
    instrument = prod.instrument.lower()
    grating = prod.grating.lower()
    target = prod.target
    name = "hlsp_ullyses_hst_{}_{}_{}_{}_cspec.fits".format(instrument, target, grating, version)
    return name

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default="/astro/ullyses/ULLYSES_DATA/",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Directory for output HLSPs")
    args = parser.parse_args()

    main(args.indir, args.outdir)
