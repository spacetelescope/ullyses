VERSION = "dr4"
CAL_VER = 1.1

CUSTOM_DIR = "/astro/ullyses/custom_cal"
DATA_DIR = "/astro/ullyses/ULLYSES_DATA"
HLSP_DIR = "/astro/ullyses/ULLYSES_HLSP"
VETTED_DIR = f"/astro/ullyses/all_vetted_data_{VERSION}"

# Some targets have periods in their name and these can break MAST ingest
# Rename them to remove periods and strip any trailing numbers after periods
RENAME = {"moa-j010321.3-720538": "moa-j010321-720538",
          "sstc2dj160830.7-382827": "sstc2dj160830-382827",
          "echa-j0844.2-7833": "echa-j0844-7833",
          "sstc2dj160000.6-422158": "sstc2dj160000-422158",
          "echa-j0843.3-7915": "echa-j0843-7915",
          "ogle-j004942.75-731717.7": "ogle-j004942-731717",
          "sstc2dj161344.1-373646": "sstc2dj161344-373646"}
