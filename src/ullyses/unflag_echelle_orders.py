"""
The lowest echelle orders of STIS E140M/1425, E230M/2707, and E230M/2415 
that fit on the STIS MAMA detectors have recently been restored to the 
RIPTAB.  Although these orders are routinely extracted, the 11th DQ bit 
(2**11 = 2048; this indicates that > 30% of the background pixels are 
rejected) is always set because the background regions were not adjusted, 
and one of the regions used to estimate the background always extends 
beyond the edge of the detector.   As a result, the additional wavelength 
coverage provided by E140M order 86, E230M/2707 order 66, and E230M/2415 
order 73 is excluded from ULLYSES HLSPs during SDQ screening.  
This is particularly undesirable in the case of E140M order 86, which 
includes the important N IV 1718 feature.
"""

from astropy.io import fits
import numpy as np

def unflag_2048(filename):
    orders = {1425: 86, 2415: 73, 2707: 66} #cewnave/order pairs
    opt_elem = fits.getval(filename, "opt_elem")
    if opt_elem not in ["E140M", "E230M"]:
        print(f"Filename {filename} is not in affected modes (E140M, E230M), skipping")
        return
    cenwave = fits.getval(filename, "cenwave")
    if cenwave not in orders:
        print(f"Filename {filename} is not in affected cenwaves {tuple(orders.keys())}, skipping")
        return
    order = orders[cenwave]
    performed = True
    with fits.open(filename, mode="update") as hdulist:
        for i in range(len(hdulist)):
            if hdulist[i].name == "SCI":
                ind = np.where(hdulist[i].data["sporder"] == order)
                if len(ind[0]) == 0:
                    performed = False
                    continue
                hdulist[i].data["DQ"][ind[0][0]] -= np.bitwise_and(hdulist[i].data["DQ"][ind[0][0]], 2048)
        # Since custom processing was performed, mark these as level0
        if performed is True:
            hdulist[0].header["HLSP_LVL"] = 0
            print(f"Removed all DQ=2048 flags from {opt_elem}/{cenwave} order={order} for {filename}")
