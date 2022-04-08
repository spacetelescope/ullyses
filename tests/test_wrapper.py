# Test the wrapper script
import os
import glob
import shutil
import urllib.request
import tarfile
import pytest

from astropy.io.fits import FITSDiff

from ullyses import wrapper
from ullyses_utils.ullyses_config import VERSION
import ullyses_utils

ULLYSES_DATA_LOCATION = 'https://ullyses.stsci.edu/files/'

COS_TARGETS = ['sk-70d16',
               'lh9-89']

STIS_TARGETS = ['lh9-34',
                'pgmw3120']

MULTIPLE_TARGETS = ['sk-67d22',
                    'sk-68d26']

class TestWrapper():

    def test_cos_data(self):
        for target in COS_TARGETS:
            self.setup_tree(target)
            self.run_wrapper('./')
            self.compare_outputs()
            self.cleanup(target)
        return

    def test_stis_data(self):
        for target in STIS_TARGETS:
            self.setup_tree(target)
            self.run_wrapper('./')
            self.compare_outputs()
            self.cleanup(target)
        return

    def test_fuse_data(self):
        pass

    def test_multiple_data_types(self):
        for target in MULTIPLE_TARGETS:
            self.setup_tree(target)
            self.run_wrapper('./')
            self.compare_outputs()
            self.cleanup(target)
        return

    def setup_tree(self, target):
        if os.path.isdir(target):
            shutil.rmtree(target)
        filename = target + '.tar.gz'
        fullurl = ULLYSES_DATA_LOCATION + filename
        print("Retrieving {}".format(fullurl))
        urllib.request.urlretrieve(fullurl, filename)
        tarball = tarfile.open(filename)
        tarball.extractall()
        os.chdir(target)
        os.mkdir('truth')
        hlsps = glob.glob('hlsp_ullyses*')
        for datafile in hlsps:
            os.rename(datafile,'truth/'+datafile)
        return

    def run_wrapper(self, indir):
        wrapper.main(indir, outdir=indir, version=VERSION)
        return

    def compare_outputs(self):
        # Outputs from current run are in ./, truth files to compare
        # with are in ./truth
        new_hlsps = glob.glob('hlsp_ullyses*')
        all_ok = True
        fitsdiff_report = ''
        keywords_to_ignore = ['DATE', 'FITS_SW']
        for new_product in new_hlsps:
            truth_file = './truth/' + new_product
            fdiff = FITSDiff(new_product, truth_file, ignore_hdus=['provenance'],
            ignore_keywords=keywords_to_ignore,
            rtol=1.0e-7)
            fitsdiff_report += fdiff.report()
            if not fdiff.identical and all_ok:
                all_ok = False
        if not all_ok:
            raise AssertionError(os.linesep + fitsdiff_report)
        return

    def cleanup(self, target):
        os.chdir('..')
        shutil.rmtree(target)
        os.remove(target+'.tar.gz')
        return
