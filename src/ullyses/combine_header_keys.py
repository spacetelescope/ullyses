from abc import ABC
import os
import pandas as pd
import numpy as np
from datetime import datetime as dt
from astropy.time import Time

import ullyses_utils

SECONDS_PER_DAY = 86400.
UTILS_DIR = ullyses_utils.__path__[0]

class KeyBlender(ABC):
    def combine_keys(self, key, method="multi", dict_key=None, constant=None):
        """
        Combine keyword values from multiple input files.
        Input:
            self (instance): Instance of either a Ullyses class or SegmentList class.
            key (str): keyword that is to be combined- this is the output keyword name
            method (str): Default=multi. Method of combining keywords. Allowed values are: 
                multi (returns the value MULTI if input files have different key vals),
                min (returns minimum of all input values),
                max (returns maximum of all input values),
                average (returns average of all input values),
                sum (returns sum of all input values),
                arr (returns numpy array of all input values),
                concat (returns pipe-separated string of all input values)
            dict_key (str): Default=None. Telescope/instrument to look up exact keyword 
                location from. If None, value is looked up on the fly.
            constant (str): Default=None. If not None, this value is returned.
    
        Returns:
            Using the supplied method, a single value is returned that distills all input
                file values into one descriptor.
        """
    
        keymap= {"HST": {"expstart": ("expstart", 1),
                         "expend": ("expend", 1),
                         "exptime": ("exptime", 1),
                         "telescop": ("telescop", 0),
                         "instrume": ("instrume", 0),
                         "detector": ("detector", 0),
                         "opt_elem": ("opt_elem", 0),
                         "filter": ("filter", 0),
                         "fgslock": ("fgslock", 0),
                         "gyromode": ("gyromode", 0),
                         "flashdur": ("flashdur", 0),
                         "flashcur": ("flashcur", 0),
                         "flashlvl": ("flashlvl", 0),
                         "flashsta": ("flashsta", 0),
                         "cenwave": ("cenwave", 0),
                         "aperture": ("aperture", 0),
                         "obsmode": ("obsmode", 0),
                         "proposid": ("proposid", 0),
                         "centrwv": ("centrwv", 0),
                         "minwave": ("minwave", 0),
                         "maxwave": ("maxwave", 0),
                         "filename": ("filename", 0),
                         "specres": ("specres", 0),
                         "comment": ("special", 0),
                         "equinox": ("equinox", 0),
                         "cal_ver": ("cal_ver", 0)},
                "WFC3": {"expstart": ("expstart", 0),
                         "expend": ("expend", 0),
                         "comment": ("special", 0),
                         "equinox": ("equinox", 0),
                         "exptime": ("exptime", 0)},
                "FUSE": {"expstart": ("obsstart", 0),
                         "expend": ("obsend", 0),
                         "exptime": ("obstime", 0),
                         "telescop": ("telescop", 0),
                         "instrume": ("instrume", 0),
                         "detector": ("detector", 0),
                         "opt_elem": ("detector", 0),
                         "cenwave": ("centrwv", 0),
                         "aperture": ("aperture", 0),
                         "obsmode": ("instmode", 0),
                         "proposid": ("prgrm_id", 0),
                         "centrwv": ("centrwv", 0),
                         "minwave": ("wavemin", 0), 
                         "maxwave": ("wavemax", 0),
                         "filename": ("filename", 0),
                         "specres": ("spec_rp", 1),
                         "comment": ("special", 0),
                         "equinox": ("equinox", 0),
                         "cal_ver": ("cf_vers", 0)},
               "LCOGT": {"expstart": ("date-obs", 1),
                         "expend": ("exptime", 1),
                         "exptime": ("exptime", 1),
                         "telescop": ("telescop", 1),
                         "instrume": ("instrume", 1),
                         "detector": ("telescop", 1),
                         "opt_elem": ("telescop", 1),
                         "proposid": ("propid", 1),
                         "filename": ("origname", 1),
                         "filter": ("filter", 1),
                         "comment": ("special", 0),
                         "cal_ver": ("pipever", 1)},
                 "VLT": {"expstart": ("mjd-obs", 1),
                         "expend": ("exptime", 1), 
                         "exptime": ("exptime", 1),
                         "telescop": ("telescop", 1),
                         "instrume": ("instrume", 1),
                         "detector": ("hierarch eso det name", 1),
                         "opt_elem": ("seqarm", 1),
                         "s_region": ("hierarch eso ada posang", 1), # calculate this from position angle
                         "equinox": ("equinox", 1),
                         "radesys": ("radecsys", 1),
                         "proposid": ("hierarch eso obs prog id", 1),
                         "centrwv": ("arm", 1), 
                         "minwave": ("arm", 1), 
                         "maxwave": ("arm", 1),
                         "cenwave": ("arm", 1), 
                         "filename": ("extname", 1),
                         "dr_num": ("dr_num", 0),
                         "dr_date": ("dr_date", 0),
                         "specres": ("special", 1),
                         "comment": ("special", 0),
                         "aperture": ("special", 0),
                         "cal_ver": ("hierarch eso pro rec1 pipe id", 1)}
                         }
   
        if key == "comment":
            dbfile = os.path.join(UTILS_DIR, "data", "calibration_metadata", "ullyses_calibration_db.csv")
            cal_db = pd.read_csv(dbfile, keep_default_na=False)
        
        if constant is not None:
            vals = [constant for x in self.first_headers]
        else:
            vals = []
            for i in range(len(self.first_headers)):
                if dict_key is None:
                    try:
                        tel = self.telescope
                    except AttributeError:
                        tel = self.primary_headers[i]["telescop"]
                else:
                    tel = dict_key

                actual_key = keymap[tel][key][0]  
                hdrno = keymap[tel][key][1]
                
                if key == "comment":
                    filename_key = keymap[tel]["filename"][0]  
                    filename_hdrno = keymap[tel]["filename"][1]
                    if filename_hdrno == 0:
                        filename = self.primary_headers[i][filename_key]
                    else:
                        filename = self.first_headers[i][filename_key]
                    db_roots = cal_db["dataset_name"].values
                    val = ""
                    for rootname in db_roots:
                        if rootname in filename:
                            val = cal_db.loc[cal_db.dataset_name == rootname]["qual_comm"].values[0]
                            break
                elif tel == "FUSE" and key == "filename":
                    val = self.primary_headers[i][actual_key]
                    val = val.replace(".fit", "_vo.fits")
                elif tel == "LCOGT" and key == "telescop":
                    telescop = self.first_headers[i]["telescop"]
                    val = f"LCOGT-{telescop}"
                elif tel == "LCOGT" and (key == "expstart" or key == "expend"):
                    dto = dt.strptime(self.first_headers[i]["date-obs"], "%Y-%m-%dT%H:%M:%S.%f")
                    t = Time(dto, format="datetime")
                    mjdstart = t.mjd
                    if key == "expstart":
                        val = mjdstart
                    elif key == "expend":
                        exptime = self.first_headers[i]["exptime"]
                        val = mjdstart + (exptime / SECONDS_PER_DAY) 
                elif tel == "VLT" and key == "expend":
                    mjdstart = self.first_headers[i]["mjd-obs"]
                    exptime = self.first_headers[i]["exptime"]
                    val = mjdstart + (exptime / SECONDS_PER_DAY)
                elif tel == "VLT" and (key == "minwave" or key == "maxwave" or key == "cenwave"):
                    if self.first_headers[i]["arm"] == "UVB": 
                        wave_vals = {"minwave": 3000, "maxwave": 5551, "cenwave": 4276}
                    else:
                        wave_vals = {"minwave": 5451, "maxwave": 10202, "cenwave": 7827}
                    val = wave_vals[key]
                elif tel == "VLT" and key == "aperture":
                    if self.first_headers[i]["arm"] == "UVB":   
                        val = self.first_headers[i]["hierarch eso ins opti3 name"] 
                    else:
                        val = self.first_headers[i]["hierarch eso ins opti4 name"]  
                elif tel == "VLT" and key == "specres":
                    if self.first_headers[i]["arm"] == "UVB": 
                        val = 6700 
                    else:
                        val = 11400
                else:
                    if hdrno == 0:
                        val = self.primary_headers[i][actual_key]
                    else:
                        val = self.first_headers[i][actual_key]
    
    #            match [tel, key]:
    #                # Handle some special cases
    #                case ["FUSE", "filename"]:
    #                    val = val.replace(".fit", "_vo.fits")
    #                case ["LCOGT", "telescop"]:
    #                    telescop = self.first_headers[i]["telescop"] 
    #                    val = f"LCOGT-{telescop}"
    #                case ["LCOGT", "expstart" | "expend" as k]:
    #                    dto = dt.strptime(self.first_headers[i]["date-obs"], "%Y-%m-%dT%H:%M:%S.%f")
    #                    t = Time(dto, format="datetime")
    #                    mjdstart = t.mjd
    #                    if k == "expstart":
    #                        val = mjdstart
    #                    if k == "expend":
    #                        exptime = self.first_headers[i]["exptime"]
    #                        val = mjdstart + (exptime / SECONDS_PER_DAY) 
    #                case ["VLT", "expend"]:
    #                    mjdstart = self.first_headers[i]["mjd-obs"]
    #                    exptime = self.first_headers[i]["exptime"]
    #                    val = mjdstart + (exptime / SECONDS_PER_DAY) 
    #                case ["VLT", "minwave" | "maxwave" | "cenwave" as k]:
    #                    if self.first_headers[i]["arm"] == "UVB":
    #                        wave_vals = {"minwave": 3000, "maxwave": 5551, "cenwave": 4276}
    #                    else:
    #                        wave_vals = {"minwave": 5451, "maxwave": 10202, "cenwave": 7827}
    #                    val = wave_vals[k]
    #                case ["VLT", "aperture"]:
    #                    if self.first_headers[i]["arm"] == "UVB":
    #                        val = self.first_headers[i]["hierarch eso ins opti3 name"]
    #                    else:
    #                        val = self.first_headers[i]["hierarch eso ins opti4 name"]
    #                case ["VLT", "specres"]:
    #                    if self.first_headers[i]["arm"] == "UVB":
    #                        val = 6700
    #                    else:
    #                        val = 11400
    #                # For normal cases
    #                case other:
    #                    actual_key = keymap[tel][key][0]
    #                    hdrno = keymap[tel][key][1]
    #                    if hdrno == 0:
    #                        val = self.primary_headers[i][actual_key]
    #                    else:
    #                        val = self.first_headers[i][actual_key]
    
                vals.append(val)
    
        # Allowable methods are min, max, average, sum, multi, arr, concat
        if method == "multi":
            keys_set = list(set(vals))
            if len(keys_set) > 1:
                return "MULTI"
            else:
                return keys_set[0]
        elif method == "min":
            return min(vals)
        elif method == "max":
            return max(vals)                                                                 
        elif method == "average": 
            return np.average(vals)
        elif method == "sum":
            return np.sum(vals)
        elif method == "arr":
            return np.array(vals)
        elif method == "concat":
            vals = [x for x in vals if x != ""]
            return " | ".join(vals)
        elif method == "comment":
            comm_vals = []
            for comm in vals:
                if comm not in comm_vals and comm != "":
                    comm_vals.append(comm)
            comm_lines = []
            for comm in comm_vals:
                starti = 0
                stopi = 0
                line_l = 0
                sp = comm.split(" ")
                for i,word in enumerate(sp):
                    line_l += len(word)
                    if i != len(sp)-1:
                        line_l += 1
                    if line_l > 72:
                        stopi = i
                        newline = " ".join(sp[starti:stopi])
                        comm_lines.append(newline)
                        line_l = 0
                        starti = i
                        if i == len(sp)-1:
                            comm_lines.append(sp[-1])
                    elif i == len(sp)-1:
                        stopi = i+1
                        newline = " ".join(sp[starti:stopi])
                        comm_lines.append(newline)

            return comm_lines
