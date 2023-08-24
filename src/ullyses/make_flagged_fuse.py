import shutil
from astropy.io import fits as pf
import glob
import os
import numpy as np
import matplotlib.pyplot as pl

from fuse_add_dq import add_dq_col

drdir = "/astro/ullyses/all_vetted_data_dr2"

# DR1 FUSE targets that require custom flagging
filestoedit = {
"AV232":         {"minwl": [1141],      "maxwl":[-1],           "dq":[1]},
"AV15":          {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"AV16":          {"minwl": [1076],      "maxwl":[1090],         "dq":[2]},
"SK-69D104":     {"minwl": [1079,992],  "maxwl":[1090,998],     "dq":[2,2]},
"AV18":          {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1] },
"LMCX-4":        {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"N11-ELS-013":   {"minwl": [0,1082],    "maxwl":[999,1094],     "dq":[2,2]},
"AV321":         {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-70D79":      {"minwl": [1136],      "maxwl":[1164],         "dq":[1]},
"AV69":          {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-67D211":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-68D135":     {"minwl": [1150],      "maxwl":[-1],           "dq":[1]},
"AV47":          {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-69D191":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"NGC346-ELS-07": {"minwl": [0,1082.5],  "maxwl":[987.5,1094.2], "dq":[2,2]},
"AV332":         {"minwl": [0,1082.5],  "maxwl":[987.6,1094.2], "dq":[2,2]},
"SK-67D101":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-70D115":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-71D46":      {"minwl": [0,1079.5],  "maxwl":[998,1094.5],   "dq":[2,2]},
"AV362":         {"minwl": [1080],      "maxwl":[1090],         "dq":[2]},
"SK-71D41":      {"minwl": [0,1079.6],  "maxwl":[992,1094.6],   "dq":[2,2]},
"SK-67D20":      {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-67D22":      {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"AV75":          {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-67D111":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"SK-67D108":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"NGC346-ELS-26": {"minwl": [0,1080],    "maxwl":[1000,1094],    "dq":[2,2]},
"SK-67D106":     {"minwl": [0,1079.5],  "maxwl":[998.5,1090.5], "dq":[2,2]},
"SK-67D168":     {"minwl": [1089.5],    "maxwl":[1094.5],       "dq":[2]},
"SK-67D107":     {"minwl": [0,1179.9],  "maxwl":[992,-1],       "dq":[2,1]},
"SK-67D105":     {"minwl": [1179.9],    "maxwl":[-1],           "dq":[1]},
"AV210":         {"minwl": [1150],      "maxwl":[1170],         "dq":[1]},
"N11-ELS-018":   {"minwl": [0,1090],    "maxwl":[990,1094.5],   "dq":[2,2]},
"SK-69D279":     {"minwl": [1080,1180], "maxwl":[1090,-1],      "dq":[2,1]}
}

# List of all FUSE DR1 targets. Commented lines are targets that had serious
# data quality issues and will not be used in DR1.
fuse_dr1 = [ 
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
 'AV18']

targstoedit = list(filestoedit.keys())

for targ in fuse_dr1:
    vofiles0 = glob.glob(os.path.join("/astro/ullyses/fuse_data", targ, "*_vo.fits"))
    vofiles = [x for x in vofiles0 if "dqscreened" not in x]
    if len(vofiles) > 1:
        print(f"more than one VO file for {targ}, exiting")
        break
    vofile = vofiles[0]
    vofilename = os.path.basename(vofile)
    outfilename = "dqscreened_"+vofilename
    outfile = os.path.join("/astro/ullyses/fuse_data", targ, outfilename)
    if targ in targstoedit:
        pars = filestoedit[targ] 
        add_dq_col(vofile, outfile, pars["minwl"], pars["maxwl"], pars["dq"], overwrite=True)
    else:
        add_dq_col(vofile, outfile, [], [], [], overwrite=True)

    shutil.copy(outfile, os.path.join(drdir, targ.lower()))

