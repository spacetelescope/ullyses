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
               'n11-els-018']

STIS_TARGETS = ['lh9-34',
                'pgmw3120']

MULTIPLE_TARGETS = ['sk-67d22',
                    'sk-68d26']

class TestWrapper():

    def test_cos_data(self):
        for target in COS_TARGETS:
            self.setup_tree(target)
            self.run_wrapper('./')
            report = self.compare_outputs()
            self.cleanup(target)
            if report is not None:
                raise AssertionError(report)
        return

    def test_stis_data(self):
        for target in STIS_TARGETS:
            self.setup_tree(target)
            self.run_wrapper('./')
            report = self.compare_outputs()
            self.cleanup(target)
            if report is not None:
                raise AssertionError(report)
        return

    def test_fuse_data(self):
        pass

    def test_multiple_data_types(self):
        for target in MULTIPLE_TARGETS:
            self.setup_tree(target)
            self.run_wrapper('./')
            report = self.compare_outputs()
            self.cleanup(target)
            if report is not None:
                raise AssertionError(report)
        return

    def setup_tree(self, target):
        if os.path.isdir(target):
            shutil.rmtree(target)
        if os.path.isdir('truth'):
            shutil.rmtree('truth')
        filename = target + '.tar.gz'
        fullurl = ULLYSES_DATA_LOCATION + filename
        print("Retrieving {}".format(fullurl))
        urllib.request.urlretrieve(fullurl, filename)
        tarball = tarfile.open(filename)
        tarball.extractall()
        os.mkdir('truth')
        os.chdir(target)
        hlsps = glob.glob('hlsp_ullyses*')
        for datafile in hlsps:
            os.rename(datafile,'../truth/'+datafile)
        return

    def run_wrapper(self, indir):
        wrapper.main(indir, outdir=indir, version=VERSION)
        return

    def compare_outputs(self):
        report = None
        # Outputs from current run are in ./, truth files to compare
        # with are in ./truth
        new_hlsps = glob.glob('hlsp_ullyses*')
        all_ok = True
        fitsdiff_report = ''
        keywords_to_ignore = ['DATE', 'FITS_SW', 'FILENAME', 'HLSP_VER']
        for new_product in new_hlsps:
            truth_filename = self.get_truth_filename(new_product)
            fdiff = FITSDiff(new_product, truth_filename, ignore_hdus=['provenance'],
            ignore_keywords=keywords_to_ignore,
            rtol=1.0e-7)
            fitsdiff_report += fdiff.report()
            if not fdiff.identical and all_ok:
                all_ok = False
        if not all_ok:
            report = os.linesep + fitsdiff_report
            return report
        return None

    def get_truth_filename(self, product):
        # Get the truth filename.  The data release might be different
        dr_position = product.find('_dr')
        if dr_position == -1:
            print('Cannot find dr version in {}'.format(product))
            return None
        stop = dr_position + 3
        filestring = '../truth/' + product[:stop] + '*'
        truth_filename_list = glob.glob(filestring)
        if len(truth_filename_list) > 1:
            print('More than 1 truth filename matches {} != 1'.format(product))
            for file in truth_filename_list:
                print(file)
            print('Returning first instance: {}'.format(truth_filename_list[0]))
        elif len(truth_filename_list) == 0:
            print('No truth files match product specification {}'.format(filestring))
            return None
        return truth_filename_list[0]
            
        return truth_filename

    def cleanup(self, target):
        os.chdir('..')
        shutil.rmtree(target)
        shutil.rmtree('truth')
        os.remove(target+'.tar.gz')
        return
