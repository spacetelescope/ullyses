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
    :return: None
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', required=True, help='Directory with corrtag data.')
    parser.add_argument('-o', '--outdir', required=True, help='Directory for the splittag and calcos output.')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('-t', '--timelist',
                       help='String list of times to split exposures over. Ex: "0, 20, 40, 60".')
    group.add_argument('-s', '--startstopincr',
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
    :rootdir: Current working directory
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


def runsplittag(inputfile, starttime, increment, endtime, timelist):
    """
    Runs the splittag command
    Needs inputs of either start, stop, and increment or timelist to run
    :param inputfile: Corrtag filename
    :param starttime: Time to start the split
    :param increment: Time increment of split
    :param endtime: Tme to end the split
    :param timelist: String of times to split over. ex: "0, 10, 20, 30, 40"
    :return: None
    """

    f1 = fits.open(inputfile)
    rootname = f1[0].header['rootname']
    root = 'split_' + rootname
    splittag.splittag(inputfile, root,
                      starttime=starttime, increment=increment, endtime=endtime,
                      time_list=timelist)


def calibratefiles(outputfolder, item):
    """
    Runs calcos on specified files
    :param outputfolder: folder to place calcos output in
    :param item: corrtag file to calibrate
    :return: None
    """

    calcos.calcos(item, outdir=outputfolder)


def parallelcal(files, outputfolder, ncores):
    """
    Parallelizes calcos
    :param files: Corrtag files to calibrate
    :param outputfolder: Folder to place calcos output in
    :param ncores: Number of cores for parallelization
    :return: None
    """

    pool = mp.Pool(processes=ncores)
    calfunc = partial(calibratefiles, outputfolder)
    pool.map(calfunc, files)


def rmlargefiles(outputdir):
    """
    Removes the FLT and COUNT files from the output directory
    These files are large in size and are rarely used
    :param outputdir: CalCOS output directory
    :return: None
    """

    fltfiles = glob.glob(os.path.join(outputdir, '*flt*'))
    for file in fltfiles:
        os.remove(file)
    countsfiles = glob.glob(os.path.join(outputdir, '*counts*'))
    for file in countsfiles:
        os.remove(file)


if __name__ == "__main__":

    # collect the command line arguments
    indir, outdir, tlist, start, end, incr, clob, ncores = parseargs()

    # get current working directory
    cwdir = os.getcwd()

    # create the output directory if it doesn't exist,
    # and remove old output if clob=True
    setupdirs(cwdir, outdir, clob)

    # collect all the corrtags to split
    corrtagfiles = glob.glob(os.path.join(cwdir, indir, '*corrtag*'))

    # g160m = [x for x in corrtagfiles if fits.getval(x, 'opt_elem') == 'G160M']
    g230l = [x for x in corrtagfiles if fits.getval(x, 'opt_elem') == 'G230L']

    # run splittag on all the corrtags
    for corrtag in g230l:
        runsplittag(corrtag, start, incr, end, tlist)

    # splittag outputs to the cwd, so move to the right output folder
    mvsplittagoutput(cwdir, outdir)

    # collect all the split corrtag files
    splitcorrtagfiles = glob.glob(os.path.join(cwdir, outdir, '*split*'))

    # filter them to not include the b segment
    splitcorrtagfiles_nob = [x for x in splitcorrtagfiles if '_b' not in x]

    # run calcos on the split corrtags in parallel
    if ncores is not None:
        parallelcal(splitcorrtagfiles_nob, os.path.join(cwdir, outdir, 'calcosout'), ncores=ncores)

    # run calcos on the split corrtags without parallelizing
    else:
        for splitfile in splitcorrtagfiles_nob:
            calibratefiles(os.path.join(cwdir, outdir, 'calcosout'), splitfile)

    # remove the FLT and COUNTS files
    rmlargefiles(os.path.join(cwdir, outdir))
