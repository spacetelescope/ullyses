import argparse
import os
import glob
import numpy as np

from astropy.io import fits

from coadd import COSSegmentList, STISSegmentList

version = 'v0.1'

'''
This wrapper goes through each target folder in the ullyses data directory and find
the data and which gratings are present. This info is then fed into coadd.py.
'''


def main(indir, outdir, version_=version):

    for root, dirs, files in os.walk(indir, topdown=False):

        print(root)
        targetname = root.split('/')[-1]
        print(f"   {targetname}")

        # collect the gratings that we will loop through
        # coadd.py will find the correct files itself,
        # but we need to know which gratings are present
        uniqmodes = []

        for myfile in glob.glob(os.path.join(root, '*_x1d.fits')):
            f1 = fits.open(myfile)
            prihdr = f1[0].header
            obsmode = (prihdr['INSTRUME'], prihdr['OPT_ELEM'])
            if obsmode not in uniqmodes:
                uniqmodes.append(obsmode)

        if not uniqmodes:
            print(f'No data to coadd for {targetname}.')
            continue

        for instrument, grating in uniqmodes:
            products = {}
            products['g130m'] = None
            products['g160m'] = None
            products['g185m'] = None
            products['e140m'] = None
            products['e230m'] = None
            products['e140h'] = None
            products['e230h'] = None
            products['cos_uv_m'] = None
            # this instantiates the class
            if instrument == 'COS':
                prod = COSSegmentList(grating, path=root)
            elif instrument == 'STIS':
                prod = STISSegmentList(grating, path=root)
            else:
                print(f'Unknown mode [{instrument}, {grating}]')

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
            products[grating] = prod

        # Create Level 3 products by abutting level 2 products
#            products['cos_uv_m'] = coadd.abut(products['g130m'], products['g160m'])
#            products['cos_m'] = coadd.abut(products['cos_m'], products['g185m'])
#            products['stis_m'] = coadd.abut(products['e140m'], products['e230m'])
#            products['stis_h'] = coadd.abut(products['e140h'], products['e230h'])


def create_output_file_name(prod):
    instrument = prod.instrument.lower()
    grating = prod.grating.lower()
    target = prod.target
    name = "hlsp_ullyses_hst_{}_{}_{}_{}_cspec.fits".format(instrument, target, grating, version_)
    return name

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default="/astro/ullyses/ULLYSES_DATA/",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Directory for output HLSPs")
    parser.add_argument("-v", "--version", default=version, 
    					help="Version number of the HLSP"
    args = parser.parse_args()

    main(args.indir, args.outdir, version_=args.version)
