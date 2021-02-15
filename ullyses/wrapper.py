import argparse
import os
import glob
import numpy as np

from astropy.io import fits

from coadd import COSSegmentList, STISSegmentList, FUSESegmentList
from coadd import abut

default_version = 'dr2'
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

        x1dfiles = glob.glob(os.path.join(root, '*_x1d.fits'))
        vofiles = glob.glob(os.path.join(root, '*_vo.fits'))

        if x1dfiles:
            for myfile in x1dfiles:
                f1 = fits.open(myfile)
                prihdr = f1[0].header
                obsmode = (prihdr['INSTRUME'], prihdr['OPT_ELEM'])
                if obsmode not in uniqmodes:
                    uniqmodes.append(obsmode)
                f1.close()
        if vofiles:
            if len(vofiles) != 1:
                print("More than 1 FUSE data file, aborting")
            else:
                for myfile in vofiles:
                    f1 = fits.open(myfile)
                    prihdr = f1[0].header
                    obsmode = ('FUSE', 'FUSE')
                    uniqmodes.append(obsmode)
                    f1.close()
              
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
        products['all_hst'] = None
        products['fuse'] = None
        products['all'] = None

        level = 2
        for instrument, grating in uniqmodes:
            # this instantiates the class
            if instrument == 'COS':
                prod = COSSegmentList(grating, path=root)
            elif instrument == 'STIS':
                prod = STISSegmentList(grating, path=root)
            elif instrument == 'FUSE':
                prod = FUSESegmentList(grating, path=root)
                products['fuse'] = prod
            else:
                print(f'Unknown mode [{instrument}, {grating}]')

            # these two calls perform the main functions
            if len(prod.members) > 0:
                prod.create_output_wavelength_grid()
                prod.coadd()
                # this writes the output file
                # If making HLSPs for a DR, put them in the official folder
                prod.target = prod.ull_targname()
                prod.targ_ra, prod.targ_dec = prod.ull_coords()
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

        if products['fuse'] is not None:
            filename = create_output_file_name(products['fuse'], version, level=level)
            products['fuse'].write(filename, clobber, level=level, version=version)

        level = 4

        if products['cos_m'] is not None and products['stis_m'] is not None:
            products['all_hst'] = abut(products['cos_m'], products['stis_m'])
        elif products['cos_m'] is not None and products['stis_h'] is not None:
            products['all_hst'] = abut(products['cos_m'], products['stis_h'])
        elif products['stis_m'] is not None and products['fuse'] is not None:
            products['all_hst'] = products['stis_m']
        elif products['E140H'] is not None and products['fuse'] is not None:
            products['all_hst'] = products['stis_h']

        if products['all_hst'] is not None and products['fuse'] is not None:
            products['all'] = abut(products['fuse'], products['all_hst'])
            filename = create_output_file_name(products['all'], version, level=level)
            filename = os.path.join(outdir, filename)
            products['all'].write(filename, clobber, level=level, version=version)
            print(f"   Wrote {filename}")
        elif products['all_hst'] is not None:
            filename = create_output_file_name(products['all_hst'], version, level=level)
            filename = os.path.join(outdir, filename)
            products['all_hst'].write(filename, clobber, level=level, version=version)
            print(f"   Wrote {filename}")


def create_output_file_name(prod, version=default_version, level=3):
    instrument = prod.instrument.lower()   # will be either cos, stis, or fuse. If abbuted can be cos-stis or cos-stis-fuse
    grating = prod.grating.lower()
    target = prod.target.lower()
    version = version.lower()
    aperture = prod.aperture.lower()

    if level == 1:
        suffix = "mspec"
        tel = 'hst'
    elif level == 3 or level == 2:
        if instrument == 'fuse':
            tel = 'fuse'
            instrument = 'fuv'   # "fuv" is the "instrument" equivalent for fuse
            grating = aperture   # the grating for fuse data is set to "fuse" to change to use aperture
            suffix = 'cspec'
        else:
            tel= 'hst'
            suffix = "cspec"
    elif level == 4:
        suffix = "sed"
        grating = "uv"
        if 'fuse' in instrument:
            tel = 'hst-fuse'
        else:
            tel = 'hst'

    # Need to add logic for uv-opt here
    name = f"hlsp_ullyses_{tel}_{instrument}_{target}_{grating}_{version}_{suffix}.fits"
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
