# Test the wrapper script
import os
import glob
import shutil
import urllib.request
import tarfile

from astropy.io.fits import FITSDiff

from .. import ullyses_coadd_abut_wrapper as wrapper
from .. import __release__


AV456 = {'name': 'av-456',
         'url': 'https://stsci.box.com/shared/static/ztdzutldcl2phgh39rx9alj4tisnhoce.gz'}
AV321 = {'name': 'av-321',
         'url': 'https://stsci.box.com/shared/static/6mn1a77yqrpv1fr4agk8yie2vc3qdsvh.gz'}
HD104237E = {'name': 'hd-104237e',
             'url': 'https://stsci.box.com/shared/static/lxme0920d2z0lswgcclgave0jpzzy239.gz'}
V_HK_ORI = {'name': 'v-hk-ori',
            'url': 'https://stsci.box.com/shared/static/dzw9gt283sz6rnekhtq88jxgcl2tq376.gz'}

RELEASE = __release__


class TestWrapper():

    def test_av456(self):
        target = AV456
        self.setup_tree(target, RELEASE)
        self.run_wrapper(target['name'])
        report = self.compare_outputs(target['name'])
        self.cleanup(target['name'])
        if report is not None:
            raise AssertionError(report)
        return

    def test_av321(self):
        target = AV321
        self.setup_tree(target, RELEASE)
        self.run_wrapper(target['name'])
        report = self.compare_outputs(target['name'])
        self.cleanup(target['name'])
        if report is not None:
            raise AssertionError(report)
        return

    def test_hd104237e(self):
        target = HD104237E
        self.setup_tree(target, RELEASE)
        self.run_wrapper(target['name'])
        report = self.compare_outputs(target['name'])
        self.cleanup(target['name'])
        if report is not None:
            raise AssertionError(report)
        return

    def test_v_hk_ori(self):
        target = V_HK_ORI
        self.setup_tree(target, RELEASE)
        self.run_wrapper(target['name'])
        report = self.compare_outputs(target['name'])
        self.cleanup(target['name'])
        if report is not None:
            raise AssertionError(report)
        return

    def setup_tree(self, target, release):
        target_name = target['name']
        url = target['url']
        if os.path.isdir(target_name):
            shutil.rmtree(target_name)
        filename = target_name + '.tar.gz'
        _ = urllib.request.urlretrieve(url, filename)
        data_tarfile = tarfile.open(filename, mode='r|gz')
        data_tarfile.extractall()
        return

    def run_wrapper(self, target):
        indir = target + '/' + RELEASE + '/input/'
        wrapper.main(indir, outdir=indir, version=RELEASE)
        return

    def compare_outputs(self, target):
        report = None
        # Outputs from current run are in ./, truth files to compare
        # with are in ./truth
        all_ok = True
        fitsdiff_report = ''
        keywords_to_ignore = ['DATE', 'FITS_SW', 'FILENAME',
                              'HLSP_VER', 'S_REGION']
        new_hlsps = glob.glob(target + '/' + RELEASE + '/input/hlsp_ullyses*')
        print(new_hlsps)
        for new_product in new_hlsps:
            truth_filename = self.get_truth_filename(target, new_product)
            fdiff = FITSDiff(new_product, truth_filename,
                             ignore_hdus=['provenance'],
                             ignore_keywords=keywords_to_ignore,
                             rtol=1.0e-7)
            fitsdiff_report += fdiff.report()
            if not fdiff.identical and all_ok:
                all_ok = False
        if not all_ok:
            report = os.linesep + fitsdiff_report
            return report
        print(fitsdiff_report)
        return None

    def get_truth_filename(self, target, product):
        # Get the truth filename.  The data release might be different
        filename = os.path.basename(product)
        truth_filename = target + '/' + RELEASE + '/truth/' + filename
        return truth_filename

    def cleanup(self, target):
        shutil.rmtree(target)
        os.remove(target+'.tar.gz')
        return
