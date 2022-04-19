import shutil
import glob
import os
import argparse
from astropy.io import fits

from ullyses_utils.parse_csv import parse_aliases
from ullyses.fuse_add_dq import add_dq_col


"""
This script is used to flag bad wavelength regions
in FUSE data. The bad regions for each target are
listed in TARGS_TO_FLAG.


Inputs:
    indir (str): Path to input directory
    outdir (str): Path to output directory
    overwrite (bool): If true, files in outdir will be overwritted. Default=False
"""


# TARGS_TO_FLAG lists FUSE targets that require custom flagging
# DQ=1 (Worm)
# DQ=2 (Poor photometric quality)
TARGS_TO_FLAG = {
    # DR2 targets below
    "AV232": {"minwl": [1141], "maxwl": [-1], "dq": [1]},
    "AV15": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "AV16": {"minwl": [1076], "maxwl": [1090], "dq": [2]},
    "SK-69D104": {"minwl": [1079, 992], "maxwl": [1090, 998], "dq": [2, 2]},
    "AV18": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "LMCX-4": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "N11-ELS-013": {"minwl": [0, 1082], "maxwl": [999, 1094], "dq": [2, 2]},
    "AV321": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-70D79": {"minwl": [1136], "maxwl": [1164], "dq": [1]},
    "AV69": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-67D211": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-68D135": {"minwl": [1150], "maxwl": [-1], "dq": [1]},
    "AV47": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-69D191": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "NGC346-ELS-07": {"minwl": [0, 1082.5], "maxwl": [987.5, 1094.2], "dq": [2, 2]},
    "AV332": {"minwl": [0, 1082.5], "maxwl": [987.6, 1094.2], "dq": [2, 2]},
    "SK-67D101": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-70D115": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-71D46": {"minwl": [0, 1079.5], "maxwl": [998, 1094.5], "dq": [2, 2]},
    "AV362": {"minwl": [1080], "maxwl": [1090], "dq": [2]},
    "SK-71D41": {"minwl": [0, 1079.6], "maxwl": [992, 1094.6], "dq": [2, 2]},
    "SK-67D20": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-67D22": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "AV75": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-67D111": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "SK-67D108": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "NGC346-ELS-26": {"minwl": [0, 1080], "maxwl": [1000, 1094], "dq": [2, 2]},
    "SK-67D106": {"minwl": [0, 1079.5], "maxwl": [998.5, 1090.5], "dq": [2, 2]},
    "SK-67D168": {"minwl": [1089.5], "maxwl": [1094.5], "dq": [2]},
    "SK-67D107": {"minwl": [0, 1179.9], "maxwl": [992, -1], "dq": [2, 1]},
    "SK-67D105": {"minwl": [1179.9], "maxwl": [-1], "dq": [1]},
    "AV210": {"minwl": [1150], "maxwl": [1170], "dq": [1]},
    "N11-ELS-018": {"minwl": [0, 1090], "maxwl": [990, 1094.5], "dq": [2, 2]},
    "SK-69D279": {"minwl": [1080, 1180], "maxwl": [1090, -1], "dq": [2, 1]},
    # DR3 targets below
    "2DFS-999": {"minwl": [0], "maxwl": [1000], "dq": [2]},
    "AV26": {"minwl": [1150], "maxwl": [-1], "dq": [1]},
    # DR4 targets below
    "AV96": {"minwl": [1147], "maxwl": [1188], "dq": [1]},
    "SK-70D32": {"minwl": [0, 1082.5], "maxwl": [992, 1087.2], "dq": [2, 2]}
}

# For completeness, this is a list of all ULLYSES targets that
# have FUSE data, but data were not used due to data quality issues.
# It may be possible to address these issues in the future.
BAD_FUSE = [
    'AV83',
    'PGMW3120',
    'LH114-7',
    'AV235',
    'AV476',
    'NGC346-MPG-368',
    'AV388',
    'AV267',
    'AV22',
    'NGC346-MPG-355',
    'BAT99-105',
    'AV287',
    'SK-69D220',
    'NGC346-MPG-342',
    'MOA-J010321.3-720538',
    'NGC2004-ELS-26',
    'BI272',
    'NGC346-MPG-435',
    'SK-68D16']

# List of all FUSE targets by DR. Commented lines are targets that had serious
# data quality issues and were not included in the DR.
FUSE_TARGS = {
    "DR2": [
        'AV207',
        'SK-71D19',
        # 'AV83',
        # 'PGMW3120',
        'AV393',
        # 'LH114-7',
        'AV321',
        'AV47',
        'AV16',
        'SK-67D168',
        'AV80',
        'SK-67D191',
        'AV377',
        'AV242',
        'AV423',
        'SK-68D135',
        'AV304',
        'AV372',
        'SK-67D105',
        # 'AV235',
        'SK-70D115',
        'LMCX-4',
        'SK-67D5',
        'SK-67D14',
        'AV446',
        'AV232',
        'N11-ELS-013',
        'AV479',
        'SK-69D50',
        # 'AV476',
        'SK-68D140',
        'AV266',
        'AV440',
        'SK-67D2',
        'SK-65D22',
        'SK-71D46',
        'SK-67D108',
        'AV229',
        'NGC346-ELS-07',
        # 'NGC346-MPG-368',
        'SK-67D101',
        'SK-67D107',
        'AV243',
        'SK-67D22',
        # 'AV388',
        'AV210',
        'SK-67D20',
        'SK-69D104',
        'VFTS72',
        'AV175',
        'AV95',
        'SK-71D41',
        'SK-71D50',
        'SK-69D191',
        'AV362',
        'SK-68D73',
        'AV69',
        'NGC346-ELS-43',
        'SK-65D47',
        'SK-69D279',
        'BI237',
        'SK191',
        'AV70',
        'SK-67D211',
        'AV6',
        'SK-68D15',
        'AV15',
        'AV187',
        'SK-68D52',
        'SK-70D79',
        'SK-67D106',
        'BI173',
        'SK-68D26',
        'BI184',
        'SK-66D35',
        'AV104',
        'SK-67D111',
        # 'AV267',
        'NGC346-ELS-26',
        'AV215',
        'AV488',
        # 'AV22',
        # 'NGC346-MPG-355',
        'SK-68D155',
        'N11-ELS-018',
        'AV456',
        # 'BAT99-105',
        'AV332',
        'AV75',
        'AV327',
        'AV18'],
    "DR3": [
        'AV490',
        'SK-67D166',
        # 'AV287',
        'SK-66D51',
        'HD38029',
        'SK-70D60',  # This data was good, there was no accompanying HST data for DR3
        'SK-69D175',
        # 'SK-69D220',
        'SK-67D167',
        '2DFS-999',
        'AV26',
        'SK-69D246',
        'AV216',
        'SK188',
        # 'NGC346-MPG-342',
        # 'MOA-J010321.3-720538',
        # 'NGC2004-ELS-26',
        'AV177',
        'SK-66D172',
        'SK190',
        'SK-67D104',
        'SK-68D129',
        # 'BI272',
        'HV5622',
        # 'NGC346-MPG-435',
        # 'SK-68D16',
        'SK-67D118'],

    "DR4": [
        "AV220",
        "AV472",
        "AV261",
        "AV264",  # This data was good, but there was no accompanying HST data in DR4
        "AV96",  # This data was good, but there was no accompanying HST data in DR4
        "HD269927C",  # This data was good, but there was no accompanying HST data in DR4
        "SK-65D55",  # This data was good, but there was no accompanying HST data in DR4
        # "SK-66D171", # This data was good, but there was no accompanying HST data in DR4
        "SK-69D43",  # This data was good, but there was no accompanying HST data in DR4
        # "SK-70D32", # This data was good, but there was no accompanying HST data in DR4
        "SK-71D21",
        "SK-71D8",  # This data was good, but there was no accompanying HST data in DR4
        "AV170",
        "AV85",  # This data was good, but there was no accompanying HST data in DR4
        "SK-67D266",
        "SK-68D41"  # This data was good, but there was no accompanying HST data in DR4
        ],

    "DR5": []
}


def flag_file(vofile, outdir, overwrite=False):
    """
    This code steps through the targets and flags
    the DQs appropriately

    :param vofile: Path to VO file
    :param outdir: Path to output directory
    :param overwrite: If true, overwrite output. Default=False
    :return: None
    """

    if "dqscreened_" in vofile:
        print(f"Input file {vofile} is already flagged, skipping")
        return vofile

    fuse_targname = fits.getval(vofile, "targname")
    aliases = parse_aliases()
    alias_mask = aliases.apply(lambda row: row.astype(str).str.fullmatch(re.escape(targ.upper())).any(), axis=1)
    if set(alias_mask) != {False}:
        ull_targname = aliases[alias_mask]["ULL_MAST_name"].values[0]
    else:
        raise KeyError(f"FUSE target {fuse_targname} not found in ULLYSES alias list")

    vofilename = os.path.basename(vofile)
    outfilename = "dqscreened_" + vofilename
    outfile = os.path.join(outdir, outfilename)
    if os.path.exists(outfile) and overwrite is True:
        os.remove(outfile)
    elif os.path.exists(outfile) and overwrite is False:
        print(f"Skipping {vofile}, output already exists and overwrite is False")
        return outfile

    targstoedit = list(TARGS_TO_FLAG.keys())
    if ull_targname in targstoedit:
        pars = TARGS_TO_FLAG[ull_targname]
        add_dq_col(vofile, outfile, pars["minwl"], pars["maxwl"], pars["dq"], overwrite=overwrite)
    else:
        print(f"No bad regions in {vofile}, but still adding DQ array")
        add_dq_col(vofile, outfile, [], [], [], overwrite=True)


def main(indir, outdir, overwrite):
    """
    Perform flagging.
    :param indir: Path to input directory
    :param outdir: Path to output directory
    :param overwrite: If true, overwrite output. Default=False
    :return: None
    """

    flag_file(indir, outdir, overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir", default=".",
                        help="Path to input directory with input NVO files to flag")
    parser.add_argument("-o", "--outdir",
                        help="Output path to place flagged NVO files")
    parser.add_argument("-c", "--overwrite", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")

    args = parser.parse_args()

    main(args.indir, args.outdir, args.overwrite)
