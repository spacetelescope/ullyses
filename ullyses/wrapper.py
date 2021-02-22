import argparse
import os
import glob
import numpy as np

from astropy.io import fits

from coadd import COSSegmentList, STISSegmentList
from coadd import abut

default_version = 'dr1'
PROD_DIR = "/astro/ullyses/ULLYSES_HLSP"

'''
This wrapper goes through each target folder in the ullyses data directory and find
the data and which gratings are present. This info is then fed into coadd.py.
'''


def main(indir, outdir, version=default_version, clobber=False):
    outdir_inplace = False
    if outdir is None:
        outdir_inplace = True
    for root, dirs, files in os.walk(indir, topdown=False):
        # Given a dir structure as follow, setting depth=2 ensure subdir/ will not be read
        # ULLYSES_DATA/
        # |___targ1/
        #     |___subdir/
        # 
        depth = 2
        if root[len(indir):].count(os.sep) >= depth:
            continue 
        
        print(root)
        dirname = root.split('/')[-1]
        print(f"   {dirname}")

        # collect the gratings that we will loop through
        # coadd.py will find the correct files itself,
        # but we need to know which gratings are present
        uniqmodes = []

        spec1d = glob.glob(os.path.join(root, '*_x1d.fits')) + glob.glob(os.path.join(root, '*_sx1.fits'))
        for myfile in spec1d:
            f1 = fits.open(myfile)
            prihdr = f1[0].header
            obsmode = (prihdr['INSTRUME'], prihdr['OPT_ELEM'])
            if obsmode not in uniqmodes:
                uniqmodes.append(obsmode)

        if not uniqmodes:
            print(f'No data to coadd for {dirname}.')
            continue

        products = {}
        products['G130M'] = None
        products['G160M'] = None
        products['G185M'] = None
        products['E140M'] = None
        products['E230M'] = None
        products['E140H'] = None
        products['E230H'] = None
        products['cos_fuv_m'] = None
        products['cos_m'] = None
        products['stis_m'] = None
        products['stis_h'] = None
        products['all'] = None

        level = 2
        for instrument, grating in uniqmodes:
            # this instantiates the class
            if instrument == 'COS':
                prod = COSSegmentList(grating, path=root)
            elif instrument == 'STIS':
                prod = STISSegmentList(grating, path=root)
            else:
                print(f'Unknown mode [{instrument}, {grating}]')
                continue

            prod.target = prod.ull_targname()
            prod.targ_ra, prod.targ_dec = prod.ull_coords()

            # these two calls perform the main functions
            if len(prod.members) > 0:
                prod.create_output_wavelength_grid()
                prod.coadd()
                # this writes the output file
                # If making HLSPs for a DR, put them in the official folder
                target = prod.target.lower()
                if outdir_inplace is True:
                    outdir = os.path.join(PROD_DIR, target, version)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                outname = create_output_file_name(prod, version, level=level)
                outname = os.path.join(outdir, outname)
                prod.write(outname, clobber, level=level, version=version)
                print(f"   Wrote {outname}")
                products[grating] = prod
            if prod.level0 is True:
                prod.create_output_wavelength_grid()
                prod.coadd()
                # this writes the output file
                # If making HLSPs for a DR, put them in the official folder
                target = prod.target.lower()
                if outdir_inplace is True:
                    outdir = os.path.join(PROD_DIR, target, version)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                outname = create_output_file_name(prod, version, level=0)
                outname = os.path.join(outdir, outname)
                prod.write(outname, clobber, level=0, version=version)
                print(f"   Wrote {outname}")
                products[grating] = prod
            else:
                print(f"No valid data for grating {grating}")
            products[grating] = prod

        # Create Level 3 products by abutting level 2 products
#            products['cos_fuv_m'] = coadd.abut(products['g130m'], products['g160m'])
#            products['cos_m'] = coadd.abut(products['cos_m'], products['g185m'])
#            products['stis_m'] = coadd.abut(products['e140m'], products['e230m'])
#            products['stis_h'] = coadd.abut(products['e140h'], products['e230h'])

        # Create Level 3 products by abutting level 2 products
        level = 3
        if products['G130M'] is not None and products['G160M'] is not None:
            products['cos_fuv_m'] = abut(products['G130M'], products['G160M'])
            filename = create_output_file_name(products['cos_fuv_m'], version, level=level)
            filename = os.path.join(outdir, filename)
            products['cos_fuv_m'].write(filename, clobber, level=level, version=version)
            print(f"   Wrote {filename}")
        elif products['G130M'] is not None:
            products['cos_fuv_m'] = products['G130M']
        elif products['G160M'] is not None:
            products['cos_fuv_m'] = products['G160M']

        if products['cos_fuv_m'] is not None and products['G185M'] is not None:
            products['cos_m'] = abut(products['cos_fuv_m'], products['G185M'])
            if products['cos_m'] is not None:
                filename = create_output_file_name(products['cos_m'], version, level=level)
                filename = os.path.join(outdir, filename)
                products['cos_m'].write(filename, clobber, level=level, version=version)
                print(f"   Wrote {filename}")
        elif products['cos_fuv_m'] is not None:
            products['cos_m'] = products['cos_fuv_m']
        elif products['G185M'] is not None:
            products['cos_m'] = products['G185M']
        
        if products['E140M'] is not None and products['E230M'] is not None:
            products['stis_m'] = abut(products['E140M'], products['E230M'])
            if products['stis_m'] is not None:
                filename = create_output_file_name(products['stis_m'], version, level=level)
                filename = os.path.join(outdir, filename)
                products['stis_m'].write(filename, clobber, level=level, version=version)
                print(f"   Wrote {filename}")
        elif products['E140M'] is not None:
            products['stis_m'] = products['E140M']
        elif products['E230M'] is not None:
            products['stis_m'] = products['E230M']
        
        if products['E140H'] is not None and products['E230H'] is not None:
            products['stis_h'] = abut(products['E140H'], products['E230H'])
            if products['stis_h'] is not None:
                filename = create_output_file_name(products['stis_h'], version, level=level)
                filename = os.path.join(outdir, filename)
                products['stis_h'].write(filename, clobber, level=level, version=version)
                print(f"   Wrote {filename}")
        elif products['E140H'] is not None:
            products['stis_h'] = products['E140H']
        elif products['E230H'] is not None:
            products['stis_h'] = products['E230H']

        level = 4
        if products['cos_m'] is not None and products['stis_m'] is not None:
            products['all'] = abut(products['cos_m'], products['stis_m'])
        elif products['cos_m'] is not None and products['stis_h'] is not None:
            products['all'] = abut(products['cos_m'], products['stis_h'])
        if products['all'] is not None:
            filename = create_output_file_name(products['all'], version, level=level)
            filename = os.path.join(outdir, filename)
            products['all'].write(filename, clobber, level=level, version=version)
            print(f"   Wrote {filename}")


def create_output_file_name(prod, version=default_version, level=3):
    instrument = prod.instrument.lower()
    grating = prod.grating.lower()
    target = prod.target.lower()
    version = version.lower()
    if level == 0:
        suffix = "spec"
    if level == 1:
        suffix = "mspec"
    elif level == 3 or level == 2:
        suffix = "cspec"
    elif level == 4:
        suffix = "sed"
        grating = "uv"
        # Need to add logic for uv-opt here
    name = f"hlsp_ullyses_hst_{instrument}_{target}_{grating}_{version}_{suffix}.fits"
    return name

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default="/astro/ullyses/all_vetted_data/",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Directory for output HLSPs")
    parser.add_argument("-v", "--version", default=default_version, 
    					help="Version number of the HLSP")
    parser.add_argument("-c", "--clobber", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")
    args = parser.parse_args()

    main(args.indir, args.outdir, args.version, args.clobber)
