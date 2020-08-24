import os
import glob
import numpy as np

from astropy.io import fits

from coadd import COSSegmentList

'''
Currently coadd.py only works for COS and only combines obs of the same grating.
This wrapper goes through each target folder in the ullyses data directory and find
the COS data and which gratings are present. This info is then feeded into coadd.py.
'''

ULLYSESDIR = '/astro/ullyses/ULLYSES_DATA/'
OUTPUTDIR = '/user/efrazer/ullyses/ullyses_dp/high_level_science_products/high_level_science_products/testcoadds/'


def main():

	for root, dirs, files in os.walk(ULLYSESDIR, topdown=False):

		print(root)
		targetname = root.split('/')[-1]
		print(targetname)

		# collect the gratings that we will loop through
		# coadd.py will find the correct files itself,
		# but we need to know which gratings are present
		modes = []

		for myfile in glob.glob(os.path.join(root, 'l*_x1d.fits')):  # only grabbing COS with the l*_x1d.fits
			f1 = fits.open(myfile)
			prihdr = f1[0].header
			modes.append(prihdr['OPT_ELEM'])

		if not modes:
			print(f'No COS data to coadd for {targetname}.')
			continue

		uniqmodes = np.unique(modes)

		for grating in uniqmodes:

			# this initiates the class
			prod = COSSegmentList(grating, path=root)

			# these two calls perform the main functions
			prod.create_output_wavelength_grid()
			prod.coadd()

			# this writes the output file
			prod.write(os.path.join(OUTPUTDIR, targetname + '_' + grating + '.fits'))


if __name__ == '__main__':

	main()
