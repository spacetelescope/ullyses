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

'''
curl -L -X GET "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast%3AHLSP%2Fullyses%2Fsk-70d16%2Fdr6%2Fhlsp_ullyses_hst_cos_sk-70d16_g130m-g160m-g185m_dr6_cspec.fits" --output "MAST_2023-06-29T15_31_07.292Z/ULLYSES/SK-70D16/hlsp_ullyses_hst_cos_sk-70d16_g130m-g160m-g185m_dr6_cspec.fits" --create-dirs

curl -L -X GET "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast%3AHLSP%2Fullyses%2Fsk-70d16%2Fdr6%2Fhlsp_ullyses_hst_cos_sk-70d16_g130m_dr6_cspec.fits" --output "MAST_2023-06-29T15_31_07.292Z/ULLYSES/SK-70D16/hlsp_ullyses_hst_cos_sk-70d16_g130m_dr6_cspec.fits" --create-dirs

curl -L -X GET "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast%3AHLSP%2Fullyses%2Fsk-70d16%2Fdr6%2Fhlsp_ullyses_hst_cos_sk-70d16_g160m_dr6_cspec.fits" --output "MAST_2023-06-29T15_31_07.292Z/ULLYSES/SK-70D16/hlsp_ullyses_hst_cos_sk-70d16_g160m_dr6_cspec.fits" --create-dirs

curl -L -X GET "https://mast.stsci.edu/api/v0.1/Download/file?uri=mast%3AHLSP%2Fullyses%2Fsk-70d16%2Fdr6%2Fhlsp_ullyses_hst_cos_sk-70d16_g185m_dr6_cspec.fits" --output "MAST_2023-06-29T15_31_07.292Z/ULLYSES/SK-70D16/hlsp_ullyses_hst_cos_sk-70d16_g185m_dr6_cspec.fits" --create-dirs

curl -L -X GET "https://mast.stsci.edu/search/hst/api/v0.1/retrieve_product?product_name=LE9R4C010%2Fle9r4c010_asn.fits" --output "MAST_2023-06-15T15_37_17.804Z/HST/LE9R4C010/le9r4c010_asn.fits" --create-dirs

curl -L -X GET "https://mast.stsci.edu/search/hst/api/v0.1/retrieve_product?product_name=LE9R4C010%2Fle9r4cx1q_x1d.fits"

'''

COS_TARGETS = {'sk-70d16': {'truth': ['hlsp_ullyses_hst_cos_sk-70d16_g130m-g160m-g185m_dr6_cspec.fits',
                                      'hlsp_ullyses_hst_cos_sk-70d16_g130m_dr6_cspec.fits',
                                      'hlsp_ullyses_hst_cos_sk-70d16_g160m_dr6_cspec.fits',
                                      'hlsp_ullyses_hst_cos_sk-70d16_g185m_dr6_cspec.fits'],
                            'input': ['le9r4cwyq_x1d.fits', 'le9r4cx1q_x1d.fits',
                                      'le9r4cxeq_x1d.fits', 'le9r4cx7q_x1d.fits',
                                      'le9r4cx4q_x1d.fits', 'le9r4cxaq_x1d.fits',
                                      'le9r4cxpq_x1d.fits', 'le9r4cxtq_x1d.fits',
                                      'le9r4cxxq_x1d.fits', 'le9r4cxzq_x1d.fits',
                                      'le9r4cxiq_x1d.fits', 'le9r4cxmq_x1d.fits',
                                      'le9r4cy2q_x1d.fits', 'le9r4cxkq_x1d.fits']},
               'n11-els-018': {'truth': [],
                               'input': []}}

STIS_TARGETS = ['lh9-34',
                'pgmw3120']

MULTIPLE_TARGETS = ['sk-67d22',
                    'sk-68d26']

RELEASE = 'dr6'

class TestWrapper():

    def test_cos_data(self):
        for target in COS_TARGETS:
            self.setup_tree(COS_TARGETS[target])
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

    def setup_tree(self, target, targetdict, release):
        if os.path.isdir(target):
            shutil.rmtree(target)
        if os.path.isdir('truth'):
            shutil.rmtree('truth')
        for ufile in targetdict['truth']:
            fullurl = self.construct_ullyses_url(target, ufile, release)
            print("Retrieving {}".format(fullurl))
            urllib.request.urlretrieve(fullurl, ufile)
#        os.mkdir('truth')
#        os.chdir(target)
#        hlsps = glob.glob('hlsp_ullyses*')
#        for datafile in hlsps:
#            os.rename(datafile,'../truth/'+datafile)
        return

    def construct_ullyses_url(self, target, filename, release):
        firstpart = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast%3AHLSP%2Fullyses%2F'
        fullurl = firstpart + target + '%2F' + release + '%2F' + filename
        return fullurl

    def construct_rawdata_url(self, filename):
        asn_id = fits.getval(filename, 'ASN_ID')
        fullurl = "https://mast.stsci.edu/search/hst/api/v0.1/retrieve_product?product_name=" + asn_id + "%2F" + filename
        return fullurl

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
