## Installation

> :warning: **This release is intended for use only in reprocessing ULLYSES Data Release 1.** :warning:
> 
> Please use the release best suited for your needs.



## ULLYSES high level science products (HLSPs)

`coadd.py` currently only works for COS and only combines observations of the same grating. It takes as inputs the grating to be combined and the path to the data. The code collects all the x1d files in the directory it's pointed to, and then filters the list to only use the files of the grating specified.
 
`wrapper.py` goes through each target in a specified folder, finds the COS data, and determines which gratings are present. This info is then fed into `coadd.py` to create the data products for each target and grating. Currently this code will collect any files in each target folder, even if the observations failed. We can probably use the tech team's csv file with the "dp acknowledgement" column to filter out good and bad data.

### Creating HLSPs 

To run the wrapper script from the command line, supply an input directory 
with COS files to combine, and an output directory to write HLSPs to.
 
    > python wrapper.py -i /path/to/data -o /output/path
    
If unspecified, these parameters will default to
`/astro/ullyses/ULLYSES_DATA` and the current directory, respectively. **NOTE**: All 
subdirectories in the input directory will be searched for COS data.

To run `coadd.py` manually:

    from coadd import COSSegmentList

    prod = COSSegmentList(grating, path='/path/to/data/')
    prod.create_output_wavelength_grid()
    prod.coadd()
    prod.write('/path/to/outputfile.fits')
    
The latest changes will still use the same methods to create HLSPs, but there are a couple of
changes in the products themselves:

 - the data are now in arrays in 1 table row, rather than having a row for each datapoint.  So when accessing
the data you need to do
```python
    f = fits.open(filename)
    table = f[1].data
    row = table[0]
    wavelength = row['wavelength']
    flux = row['flux']
```
