import os
import glob
import numpy as np
import calcos
from costools import splittag
import argparse
from astropy.io import fits
from functools import partial
import multiprocessing as mp

"""
This script uses the costools.splittag() and calcos.calcos() modules to create custom
timeseries data products for ULLYSES T-Tauri monitoring stars. Splittag divides the input
corrtags into several split-corrtags based on user-specified time intervals. These split-
corrtag files are then re-calibrated using calcos to output new x1d files. At the end, FLT and
COUNTS images (which are also output from calcos) are removed due to their large file size.

To run from the command line, users must specify:

(1) "-i" or "--indir": the input directory where the un-split corrtags reside
(2) "-o" or "--outdir": the output directory for the splittag and calcos output files
(3) the time interval(s) over which to split the corrtag, which must be formatted in 1 of 3 ways:
    (a) "-t" or "--timelist"" a string list of times in seconds to split exposures 
        over, such as "0, 20, 40, 60"
    (b) "-s" or "--startstopincr": a string list of the start and stop time and increment
        (each in seconds) to split exposures over, such as "0, 250, 50"
    (c) "-r" or "-increment": the increment of the time split in seconds, where the start 
        time=0 seconds and stop time=1000 seconds by default

Optional command line arguments are:

(1) "-c" or "--clobber": if used, the script will delete existing products in the outdir folder
(2) "-p" or "--prefix": the prefix to prepend to the split corrtags
(3) "-n" or "--ncores": if specified, the script will parallelize calcos using ncores

"""


def parseargs():
    """
    Parse command line options
    :return: None
    """

    parser = argparse.ArgumentParser()

    # i/o
    parser.add_argument('-i', '--indir', required=True, help='Directory with corrtag data.')
    parser.add_argument('-o', '--outdir', required=True, help='Directory for the splittag and calcos output.')
    parser.add_argument('-c', '--clobber', default=False, required=False,
                        action='store_true',
                        help='If True, delete existing products in outdir.')
    parser.add_argument('-p', '--prefix', default='split_',
                        help='Prefix to prepend to split corrtags')

    # splittag options
    parser.add_argument('-t', '--timelist', required=False,
                        help='String list of times to split exposures over. Ex: "0, 20, 40, 60".')
    parser.add_argument('-s', '--startstopincr', required=False,
                        help='String list of start time, stop time, and increment to split '
                             'exposures over. Ex: "0, 250, 50".')
    parser.add_argument('-r', '--increment', required=False, type=np.int,
                        help='Increment of time split. Start time=0 sec and stop time=1000 sec.')

    # parallelization
    parser.add_argument('-n', '--ncores', required=False, type=np.int,
                        help='If specified, code will parallelize calcos with n cores.')

    args = parser.parse_args()

    splitargs = [args.timelist, args.startstopincr, args.increment]
    if splitargs.count(None) != 2:
        raise IOError('Must specify one of: timelist, start/stop/increment, or split increment for splittag.')

    if args.startstopincr is not None:

        if "," not in args.startstopincr:
            raise IOError('Start/Stop/Increment must be comma-separated string, Ex: "0, 250, 50".')

        startstopincr = args.startstopincr.split(",")
        starttime = startstopincr[0]
        stoptime = startstopincr[1]
        increment = startstopincr[2]

    elif args.increment is not None:

        starttime = 0
        stoptime = 1000
        increment = args.increment

    else:
        starttime = None
        stoptime = None
        increment = None

    return args.indir.strip(), args.outdir.strip(),\
           args.timelist, starttime, stoptime, increment,\
           args.clobber, args.ncores, args.prefix

def clobberfiles(outputfolder):
    """
    Remove old output files
    :param outputfolder: Output directory to remove files from
    :return: None
    """

    rmfiles = glob.glob(os.path.join(outputfolder, 'calcosout', '*.fits')) +\
              glob.glob(os.path.join(outputfolder, '*.fits'))
    if len(rmfiles):
        print('Clobbering old output from {}'.format(outputfolder))
        for file in rmfiles:
            os.remove(file)


def setupdirs(rootdir, outputfolder, clobber):
    """
    Checks if output directory exists and if not, creates it
    Clobbers old files in output directory if specified
    :param rootdir: Current working directory
    :param outputfolder: User specified output directory
    :param clobber: Boolean switch to clean directory or not
    :return: None
    """

    if not os.path.isdir(outputfolder):
        os.mkdir(os.path.join(rootdir, outputfolder))

    if clobber:
        clobberfiles(os.path.join(rootdir, outputfolder))


def mvsplittagoutput(rootdir, outputfolder):
    """
    The splittag output is by default placed into the current working directory
    This moves the output into the specified output directory
    :param rootdir: Current working directory
    :param outputfolder: Output directory
    :return: None
    """

    outfiles = glob.glob(os.path.join(rootdir, '*split*.fits'))

    for myfile in outfiles:
        os.replace(myfile, myfile.replace(rootdir, os.path.join(rootdir, outputfolder)))


def runsplittag(inputfile, starttime, increment, endtime, timelist, outdir, 
                prefix="split_"):
    """
    Runs the splittag command
    Needs inputs of either start, stop, and increment or timelist to run
    :param inputfile: Corrtag filename
    :param starttime: Time to start the split
    :param increment: Time increment of split
    :param endtime: Tme to end the split
    :param timelist: String of times to split over. ex: "0, 10, 20, 30, 40"
    :param prefix: String prefix to prepend to split corrtags
    :return: None
    """

    f1 = fits.open(inputfile)
    rootname = f1[0].header['rootname']
    root = prefix + rootname
    path_root = os.path.join(outdir, root)
    splittag.splittag(inputfile, path_root,
                      starttime=starttime, increment=increment, endtime=endtime,
                      time_list=timelist)


def calibratefiles(outputfolder, item):
    """
    Runs calcos on specified files
    :param outputfolder: folder to place calcos output in
    :param item: corrtag file to calibrate
    :return: None
    """

    calcos.calcos(item, outdir=outputfolder, verbosity=0)


def parallelcal(files, outputfolder, ncores):
    """
    Parallelizes calcos
    :param files: Corrtag files to calibrate
    :param outputfolder: Folder to place calcos output in
    :param ncores: Number of cores for parallelization
    :return: None
    """

    with mp.Pool(ncores) as pool:
        calfunc = partial(calibratefiles, outputfolder)
        pool.map(calfunc, files)


def rmlargefiles(outputfolder):
    """
    Removes the FLT and COUNT files from the output directory
    These files are large in size and are rarely used
    :param outputdir: CalCOS output directory
    :return: None
    """

    fltfiles = glob.glob(os.path.join(outputfolder, '*flt*'))
    for file in fltfiles:
        os.remove(file)
    countsfiles = glob.glob(os.path.join(outputfolder, '*counts*'))
    for file in countsfiles:
        os.remove(file)


def main(indir, outdir, tlist=None, start=0, end=1000, incr=30, 
         clob=False, numcores=1, prefix="split_"):
    # get current working directory
    cwdir = os.getcwd()

    # create the output directory if it doesn't exist,
    # and remove old output if clob=True
    setupdirs(cwdir, outdir, clob)

    # collect all the corrtags to split
    corrtagfiles = glob.glob(os.path.join(cwdir, indir, '*corrtag*'))

    # run splittag on all the corrtags
    for corrtag in corrtagfiles:
        runsplittag(corrtag, start, incr, end, tlist, outdir, prefix)

    # splittag outputs to the cwd, so move to the right output folder
    #mvsplittagoutput(cwdir, outdir)

    # collect all the split corrtag files
    splitcorrtagfiles = glob.glob(os.path.join(cwdir, outdir, '*split*'))

    # filter them to not include the b segment
    splitcorrtagfiles_nob = [x for x in splitcorrtagfiles if '_b' not in x]

    # Avoid parallelization collisions by creating outdir a priori
    if not os.path.exists(os.path.join(outdir, 'calcosout')):
        os.makedirs(os.path.join(outdir, 'calcosout'))

    # run calcos on the split corrtags in parallel
    if numcores is not None:
        parallelcal(splitcorrtagfiles_nob, os.path.join(cwdir, outdir, 'calcosout'), ncores=numcores)

    # run calcos on the split corrtags without parallelizing
    else:
        for splitfile in splitcorrtagfiles_nob:
            calibratefiles(os.path.join(cwdir, outdir, 'calcosout'), splitfile)

    # remove the FLT and COUNTS files
    rmlargefiles(os.path.join(cwdir, outdir, 'calcosout'))

    print("\n", "~"*60, "\n", 
          f"Final calibrated products written to {os.path.join(outdir, 'calcosout')}", 
          "\n","~"*60)

if __name__ == "__main__":

    # collect the command line arguments
    indir, outdir, tlist, start, end, incr, clob, numcores, prefix = parseargs()

    # Run main function
    main(indir, outdir, tlist, start, end, incr, clob, numcores, prefix)

