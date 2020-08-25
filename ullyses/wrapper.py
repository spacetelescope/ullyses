import argparse
import os
import glob
import numpy as np

from astropy.io import fits

from coadd import COSSegmentList

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
            prod.create_output_wavelength_grid()
            prod.coadd()

            # this writes the output file
            if not os.path.exists(outdir):
                os.mkdir(outdir)
            outname = os.path.join(outdir, targetname + '_' + grating + '.fits')
            prod.write(outname)
            print(f"   Wrote {outname}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default="/astro/ullyses/ULLYSES_DATA/",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Directory for output HLSPs")
    args = parser.parse_args()

    main(args.indir, args.outdir)
