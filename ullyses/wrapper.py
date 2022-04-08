import pandas as pd
from collections import defaultdict
import argparse
import os
import glob
import numpy as np

from astropy.io import fits

from ullyses.coadd import COSSegmentList, STISSegmentList, FUSESegmentList, CCDSegmentList
from ullyses.coadd import abut
from ullyses_utils.ullyses_config import RENAME, VERSION

'''
This wrapper goes through each target folder in the ullyses data directory and find
the data and which gratings are present. This info is then fed into coadd.py.
'''

def main(indir, outdir, version=VERSION, clobber=False):
    outdir_inplace = False
    if outdir is None:
        HLSP_DIR = os.getenv('HLSP_DIR')
        if HLSP_DIR is None:
            print("Environment variable HLSP_DIR must be defined if outdir is not specified")
            raise RuntimeError("Please set HLSP_DIR and restart")
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
        vofiles = glob.glob(os.path.join(root, '*_vo.fits'))
        for myfile in spec1d:
            f1 = fits.open(myfile)
            prihdr = f1[0].header
            obsmode = (prihdr['INSTRUME'], prihdr['OPT_ELEM'], prihdr['DETECTOR'])
            if obsmode not in uniqmodes:
                uniqmodes.append(obsmode)
            f1.close()

        if vofiles:
            if len(vofiles) != 1:
                print("More than 1 FUSE data file, aborting")
            else:
                obsmode = ('FUSE', 'FUSE', 'FUSE')
                uniqmodes.append(obsmode)

        if not uniqmodes:
            print(f'No data to coadd for {dirname}.')
            continue

        # Create dictionary of all products, with each set to None by default
        products = defaultdict(lambda: None)

        level = 2
        for instrument, grating, detector in uniqmodes:
            # this instantiates the class
            if instrument == 'COS':
                prod = COSSegmentList(grating, path=root)
            elif instrument == 'STIS':
                if detector == 'CCD':
                    prod = CCDSegmentList(grating, path=root)
                else:
                    prod = STISSegmentList(grating, path=root)
            elif instrument == 'FUSE':
                prod = FUSESegmentList(grating, path=root)
                products[f'{instrument}/{grating}'] = prod
            else:
                print(f'Unknown mode [{instrument}, {grating}, {detector}]')
                continue

            prod.target = prod.ull_targname()
            prod.targ_ra, prod.targ_dec = prod.ull_coords()

            # these two calls perform the main functions
            if len(prod.members) > 0:
                prod.create_output_wavelength_grid()
                prod.coadd()
                # this writes the output file
                # If making HLSPs for a DR, put them in the official folder
                prod.target = prod.ull_targname()
                prod.targ_ra, prod.targ_dec = prod.ull_coords()
                target = prod.target.lower()
                if "." in target:
                    assert target in RENAME, f"Renaming scheme not known for {targ}"
                    dir_target = RENAME[target]
                else:
                    dir_target = target
                if outdir_inplace is True:
                    outdir = os.path.join(HLSP_DIR, dir_target, version)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                if instrument != 'FUSE': # FUSE data is written as level 3 product below
                    outname = create_output_file_name(prod, version, level=level)
                    outname = os.path.join(outdir, outname)
                    prod.write(outname, clobber, level=level, version=version)
                    print(f"   Wrote {outname}")
                products[f'{instrument}/{grating}'] = prod
            else:
                print(f"No valid data for grating {grating}")
            if prod.level0 is True:
                prod.create_output_wavelength_grid()
                prod.coadd()
                # this writes the output file
                # If making HLSPs for a DR, put them in the official folder
                target = prod.target.lower()
                if "." in target:
                    assert target in RENAME, f"Renaming scheme not known for {targ}"
                    dir_target = RENAME[target]
                else:
                    dir_target = target
                if outdir_inplace is True:
                    outdir = os.path.join(HLSP_DIR, dir_target, version)
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                outname = create_output_file_name(prod, version, level=0)
                outname = os.path.join(outdir, outname)
                prod.write(outname, clobber, level=0, version=version)
                print(f"   Wrote {outname}")
                products[f'{instrument}/{grating}'] = prod
            products[f'{instrument}/{grating}'] = prod


        # Create level 3 products- abutted spectra for gratings of the same
        # resolution for each instrument.
        level = 3
        lvl3_modes = {"cos_fuv_m": ["COS/G130M", "COS/G160M", "COS/G185M"],
                      "stis_m": ["STIS/E140M", "STIS/E230M"],
                      "stis_h": ["STIS/E140H", "STIS/E230H"],
                      "stis_l": ["STIS/G230L", "STIS/G430L", "STIS/G750L"]}
        for outprod,modes in lvl3_modes.items():
            abutted = None
            dowrite = False
            for mode in modes:
                if products[mode] is not None:
                    if abutted is None:
                        abutted = products[mode]
                    else:
                        abutted = abut(abutted, products[mode])
                        dowrite = True
            if dowrite is True:
                filename = create_output_file_name(abutted, version, level=level)
                filename = os.path.join(outdir, filename)
                abutted.write(filename, clobber, level=level, version=version)
                print(f"   Wrote {filename}")
        # Manually write out a FUSE level3 product.
        if products['FUSE/FUSE'] is not None:
            filename = create_output_file_name(products['FUSE/FUSE'], version, level=level)
            filename = os.path.join(outdir, filename)
            products['FUSE/FUSE'].write(filename, clobber, level=level, version=version)


        # Determine which gratings should contribute to the final level 4 SED HLSP.
        # Starting with the bluest product and working redward, find which products,
        # if any, overlap with the bluer product. If more than one overlaps, use
        # the one that extends further. If none overlap, still abut them- there
        # will just be a region of flux=0 in between.
        level = 4
        gratings = []
        minwls = []
        maxwls = []
        ins = []
        for instrument, grating, detector in uniqmodes:
            ins.append(instrument)
            gratings.append(grating)
            minwls.append(products[instrument+"/"+grating].first_good_wavelength)
            maxwls.append(products[instrument+"/"+grating].last_good_wavelength)
        # Only go through this exercise if there is data for more than one instrument
        if len(set(ins)) != 1:
            df = pd.DataFrame({"gratings": gratings, "ins": ins, "minwls": minwls, "maxwls": maxwls})
            used = pd.DataFrame()
            # Start with the bluest product, and remove rows from the dataframe
            # until no rows remain. The only exception is if the bluest product is
            # STIS/echelle *and* G130M+G160M combo exists. Then use G130M+G160M as bluest
            # and ignore STIS/echelle
            lowind = df["minwls"].idxmin()
            if df.loc[lowind, "gratings"] in ["E140M", "E140H"]:
                if "G130M" in gratings and "G160M" in gratings:
                    g130mind = df.loc[df["gratings"] == "G130M"].index.values
                    used = used.append(df.loc[g130mind])
                    shortestwl = df.loc[g130mind[0], "minwls"]
                    df = df.drop(index=g130mind)
                    g160mind = df.loc[df["gratings"] == "G160M"].index.values
                    used = used.append(df.loc[g160mind])
                    maxwl = df.loc[g160mind[0], "maxwls"]
                    df = df.drop(index=g160mind)
                    df = df.drop(index=lowind)
                else:
                    shortestwl = df.loc[lowind, "minwls"]
                    used = used.append(df.loc[lowind])
                    maxwl = df.loc[lowind, "maxwls"]
                    df = df.drop(lowind)
            else:
                shortestwl = df.loc[lowind, "minwls"]
                used = used.append(df.loc[lowind])
                maxwl = df.loc[lowind, "maxwls"]
                df = df.drop(lowind)
            while len(df) > 0:
                lowind = df.loc[(df["minwls"] < maxwl) & (df["maxwls"] > maxwl)].index.values
                # If G130M and G160M both exist for a given target, *always*
                # abut them together regardless of other available gratings.
                # This captures the case where there is FUSE bluer than COS/FUV.
                if "G130M" in used.gratings.values and "G160M" in gratings and "G160M" not in used.gratings.values:
                    lowind = df.loc[df["gratings"] == "G160M"].index.values
                    maxwl = df.loc[lowind[0], "maxwls"]
                    used = used.append(df.loc[lowind])
                    df = df.drop(index=lowind)
                # Handle case where more than one grating overlaps with bluer data.
                elif len(lowind) > 1:
                    df2 = df.loc[lowind]
                    ranges = df2.maxwls - df2.minwls
                    biggest = ranges.idxmax()
                    match_grating = df2.loc[biggest, "gratings"]
                    match_ind = df.loc[df["gratings"] == match_grating].index.values
                    used = used.append(df.loc[match_ind])
                    maxwl = df.loc[match_ind, "maxwls"].values[0]
                    df = df.drop(index=lowind)
                # If none overlap, abut with the next closest product.
                elif len(lowind) == 0:
                    lowind = df["minwls"].idxmin()
                    used = used.append(df.loc[lowind])
                    maxwl = df.loc[lowind, "maxwls"]
                    df = df.drop(lowind)
                # This is the easy case- only one mode overlaps with the bluer data.
                else:
                    maxwl = df.loc[lowind[0], "maxwls"]
                    used = used.append(df.loc[lowind])
                    df = df.drop(index=lowind)
                # Check every time if there are any modes that overlap completely
                # with what has been abutted so far.
                badinds = df.loc[(df["minwls"] > shortestwl) & (df["maxwls"] < maxwl)].index.values
                if len(badinds) > 0:
                    df = df.drop(index=badinds)
            # If more than one instrument was selected for abutting,
            # create level 4 product.
            if len(set(used["ins"].values)) > 1:
                abut_gr = used.iloc[0]["ins"] + "/" + used.iloc[0]["gratings"]
                abutted = products[abut_gr]
                for i in range(1, len(used)):
                    abut_gr = used.iloc[i]["ins"] + "/" + used.iloc[i]["gratings"]
                    abutted = abut(abutted, products[abut_gr])
                filename = create_output_file_name(abutted, version, level=level)
                filename = os.path.join(outdir, filename)
                abutted.write(filename, clobber, level=level, version=version)
                print(f"   Wrote {filename}")


def create_output_file_name(prod, version=VERSION, level=3):
    instrument = prod.instrument.lower()   # will be either cos, stis, or fuse. If abbuted can be cos-stis or cos-stis-fuse
    grating = prod.grating.lower()
    target = prod.target.lower()
    version = version.lower()
    aperture = prod.aperture.lower()

    # Target names can't have a period in them or it breaks MAST
    if "." in target:
        assert target in RENAME, f"Renaming scheme not known for {targ}"
        target = RENAME[target]

    if level == 0:
        tel = 'hst'
        suffix = "spec"
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
        suffix = "preview-spec"
        if "G430L" in prod.grating or "G750L" in prod.grating:
            grating = "uv-opt"
        else:
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
    parser.add_argument("-i", "--indir",
                        default=f"/astro/ullyses/all_vetted_data_{VERSION}",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Directory for output HLSPs")
    parser.add_argument("-v", "--version", default=VERSION,
                        help="Version number of the HLSP")
    parser.add_argument("-c", "--clobber", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")
    args = parser.parse_args()

    main(args.indir, args.outdir, args.version, args.clobber)
