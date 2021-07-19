import os
import calcos
from costools import splittag
import glob
import argparse
from astropy.io import fits
from functools import partial
import multiprocessing as mp

"""
This code will:
(1) use costools.splittag to divide corrtag files for the variable stars 
    into multiple corrtags based on a specified time interval
(2) run calcos on the output corrtags to re-extract them into x1ds (w/ the option to parallelize)
(3) remove unnecessary file outputs such as counts images to maintain a tidy directory

"""


def parseargs():

    """
    Parse command line options
    :return:
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', required=True, help='Directory with corrtag data.')
    parser.add_argument('-o', '--outdir', required=True, help='Directory for the splittag and calcos output.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-t', '--timelist', required=True,
                       help='String list of times to split exposures over. Ex: "0, 20, 40, 60".')
    group.add_argument('-s', '--startstopincr', required=True,
                       help='String list of start time, stop time, and increment to split '
                            'exposures over. Ex: "0, 250, 50".')

    parser.add_argument('-n', '--ncores', default=5, required=False,
                        help='If specified, code will parallelize calcos with n cores.')
    parser.add_argument('-c', '--clobber', default=False, required=False,
                        action='store_true',
                        help='If True, overwrite existing products.')

    args = parser.parse_args()

    startstopincr = args.startstopincr.split(",")
    starttime = startstopincr[0]
    stoptime = startstopincr[1]
    increment = startstopincr[2]

    return args.indir, args.outdir, args.timelist, starttime, stoptime, increment, args.clobber, args.ncores


def clobberfiles(outputdir):

    """
    Remove old output files
    :param outputdir: output directory to remove files from
    """

    rmfiles = glob.glob(os.path.join(outputdir, '*'))
    if len(rmfiles):
        print('Clobbering old output from {}'.format(outputdir))
        for file in rmfiles:
            os.remove(file)


def setupdirs(rootdir, outputdir, clobber):

    """
    Checks if output directory exists and if not, creates it
    Clobbers old files in output directory if specified
    :param outputdir: user specified output dir
    :param clobber: Boolean switch to clean dir or not
    """

    if not os.path.isdir(outputdir):
        os.mkdir(os.path.join(rootdir, outputdir))

    if clobber:
        clobberfiles(os.path.join(rootdir, outputdir))


def mvsplittagoutput(cwdir, outdir):

    outfiles = glob.glob(os.path.join(cwdir, '*split*'))

    for myfile in outfiles:

        os.replace(myfile, os.path.join())



def runsplittag(inputfile, starttime, increment, endtime, timelist):

    """
    Runs the splittag command. Needs inputs of either start, stop, and increment or timelist to run.
    :param inputfile: corrtag filename
    :param starttime: the time to start the split
    :param increment: time increment of split
    :param endtime: the time to end the split
    :param timelist: string of times to split over. ex: "0, 10, 20, 30, 40"

    :return: None
    """

    f1 = fits.open(inputfile)
    rootname = f1[0].header['rootname']
    root = 'split_' + rootname
    splittag.splittag(inputfile, root,
                      starttime=starttime, increment=increment, endtime=endtime,
                      time_list=timelist)

    mvsplittagoutput()

def calibrate_files(outputfolder, item):
    """
    Runs calcos on specified files
    :param outputfolder:
    :param item:
    :return:
    """

    calcos.calcos(item, outdir=outputfolder)


def parallel_cal(files, outputfolder, ncores):
    """
    Parallelizes calcos
    :param files:
    :param outputfolder:
    :param ncores:
    :return:
    """

    pool = mp.Pool(processes=ncores)
    calfunc = partial(calibrate_files, outputfolder)
    pool.map(calfunc, files)


def rmlargefiles(outputdir):

    fltfiles = glob.glob(os.path.join(outputdir, '*flt*'))
    for file in fltfiles:
        os.remove(file)
    countsfiles = glob.glob(os.path.join(outputdir, '*counts*'))
    for file in countsfiles:
        os.remove(file)


if __name__ == "__main__":

    indir, outdir, tlist, start, end, incr, clob, ncores = parseargs()

    cwdir = os.getcwd()

    setupdirs(cwdir, outdir, clob)

    corrtagfiles = glob.glob(os.path.join(cwdir, indir, '*corrtag*'))

    for corrtag in corrtagfiles:
        runsplittag(corrtag, start, incr, end, tlist)

    splitcorrtagfiles = glob.glob(os.path.join(cwdir, outdir, '*split*'))

    if ncores is not None:

        parallel_cal(splitcorrtagfiles, os.path.join(cwdir, outdir), ncores=ncores)

    else:

        calibrate_files(os.path.join(cwdir, outdir), splitcorrtagfiles)

    rmlargefiles(os.path.join(cwdir, outdir))
