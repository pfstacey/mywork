import sys,string,os,time,pickle,pdb,glob,astropy,os,types
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
# depreciated. Use simple_norm
#from astropy.visualization import scale_image
from astropy import units as u
import robust as rb
import djs_phot_mb as djs
from select import select
from length import *
from scipy.interpolate import interp1d as interp
import scipy as sp
import constants as c
import matplotlib.patheffects as PathEffects
import scipy.optimize as opt
from photutils import CircularAperture, SkyCircularAperture, CircularAnnulus, aperture_photometry
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#from ccdproc import cosmicray_lacosmic as la

# aperture photometry from photutils
#follow notebook

# from photutils import SkyCircularAperture, CircularAperture, CircularAnnulus
# Overplot apertures and sky annulus
# positions = [(90.73, 59.43), (73.63, 139.41), (43.62, 61.63)]
# apertures = CircularAperture(positions, r=3)
# bkg_apertures = CircularAnnulus(positions, r_in=10., r_out=15.)
# from astropy.visualization import scale_image
# plt.imshow(scale_image(data, scale='sqrt', percent=98.),origin='lower')
# apertures.plot(color='blue')
# bkg_apertures.plot(color='cyan', hatch='//', alpha=0.8)

#----------------------------------------------------------------------
# Thacher Observatory Photometry Reduction Pipeline:
#
# NOTES:
#
#    Uses the python implementation of the DJS photometry routines
#    written originally by D. Schlegel in IDL and adapted to python by
#    M. Bottom.
#
# DEPENDENCIES:
#
#    Must have astrometric routines of astrometry.net to obtain
#    astrometric solutions for data files.
#
#
# TO DO:
#
#  - Make re-centering a separate function
#  - SNR vector in "info" directory of batch_phot is not working...
#  - In "lightcurve" handle choice of reference stars more intelligently
#  - Fix inputs for lightcurve to reflect new flexibility of the routine
#  - Add aperture category to info dictionary in batch_phot
#  - Consolidate djs routines into this module
#  - Update djs routines to fit 1D and 2D Gaussian to star image.
#  - Replace use of fitgaussian with example at end of code
#
#
#  3/10/14 - Created by Jonathan Swift and Michael Bottom.
#  3/12/14 - Output directory option added for routines with output.
# (jswift) - Added readnoise calculation option in master_bias
#          - Added "done_in" time keeping routine.
#          - Fixed "choose_refs" so that it properly disconnects
#            from plotting device and returns RAs and Decs of all
#            sources
#          - Updated output from routines to minimize recalculations
#  3/17/14 - Implemented flat field option in do_phot and batch_phot
# (jswift) - Option to use specified aperture in all images in dataset
#            instead of SNR optimization.
#  6/25/14 - Added options to do_astrometry to allow for finer tuned
# (jswift)   solutions.
#          - Added option to specify skyradii for do_phot and
#            batch_phot.
#          - New annotation for bias, readnoise and dark images.
#          - Fixed subplot bug in "lightcurve" to allow for large
#            number of reference stars.
#  6/13/16 - Updated output of optimal_aperture
# (jswift)
# 12/15/16 - Cleaned up import section of module
#          - Use tqdm instead of "done_in" to monitor progress in read
#            noise calculation
#  3/24/17 - Put in check for astrometry.net installation
# (jswift) - Changed default lat/lon/elevation to Thacher Observatory
#          - Added get_info_path to check for the data diretories
#            containing the Pickles stellar templates and transmission
#            curves for throughput measurements.
#  4/21/17 - Took out dependency on astropysics. Now using astropy
#            SkyCoord.
#          - Also took out dependency on fitgaussian. Now using more
#            transparent Gaussian2D function that is internal.
#
#----------------------------------------------------------------------

# For interactive plotting
plt.ion()

import ephem

obs = ephem.Observer()
obs.lat = '34 28 00.5'
obs.lon = '-119 10 38.5'
obs.elevation = 494.7

#----------------------------------------------------------------------#
# do_astrometry:                                                       #
#----------------------------------------------------------------------#

def do_astrometry(files,clobber=False,pixlo=0.1,pixhi=1.5,ra=None,dec=None,object=None,field=0.5,
                  numstars=100,downsample=4):

    """
    Overview:
    ---------
    Input list of FITS files and return solved FITS files to same
    directory with suffix "_solved.fits"

    Requirements:
    -------------
    Astrometry.net routines with calibration tiles installed locally
    will need to manually change astrometrydotnet_dir in routine if
    not installed in default location

    Calling sequence:
    -----------------
    do_astrometry(files,clobber=True)

    """

    # Default directory for astrometry.net
    astrometrydotnet_dir = "/usr/local/astrometry"
    if not os.path.isdir(astrometrydotnet_dir):
        print astrometrydotnet_dir+' path does not exist'
        return

    for file in files:
        # Allow for graceful opt out...
        print("Press any key to quit, continuing in 1 second...")
        timeout=1
        rlist, wlist, xlist = select([sys.stdin], [], [], timeout)
        if rlist:
            break

        # Get image and header
        data, header = fits.getdata(file, 0, header=True)

        # Get telescope RA and Dec as starting point

        # Test if RA and Dec is in header
        guess = False

        if 'OBJCTRA' in header and 'OBJCTDEC' in header:
            rastr  = header['OBJCTRA']
            decstr = header['OBJCTDEC']
            coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
            RAdeg  = coords.ra.deg
            DECdeg = coords.dec.deg
            guess = True

        if object != None:
            targ = SkyCoord.from_name(object)
            RAdeg = targ.ra.degree
            DECdeg = targ.dec.degree
            guess = True

        if ra != None and dec != None:
            RAdeg = ra
            DECdeg = dec
            guess = True

        # Do some string handlings
        fname = file.split('/')[-1]
        outdir = file.split(fname)[0]+'astrometry'
        datadir = file.split(fname)[0]
        ffinal = fname.split('.')[0]+'_solved.fits'

	# Don't redo astrometry unless clobber keyword set
        if len(glob.glob(datadir+ffinal)) == 1 and not clobber:
            print("Astrometry solution for "+fname+" already exists!")
            print("Skipping...")
        else:

        # Construct the command string
            if guess:
                command=string.join(
                    [astrometrydotnet_dir+"/bin/solve-field",
                     file.rstrip(),
                     "--scale-units arcsecperpix --scale-low "+str(pixlo)+" --scale-high "+str(pixhi),
                     "--ra ",str(RAdeg)," --dec ", str(DECdeg),
                     "--radius "+str(field),
                     "--downsample "+str(downsample),
                     "--no-plots",#" --no-fits2fits ",
                     "--skip-solved",
                     "--objs "+str(numstars),
                     "--odds-to-tune-up 1e4",
                     "--no-tweak",
                     "--dir",outdir,"--overwrite"])

            else:
                command=string.join(
                    [astrometrydotnet_dir+"/bin/solve-field",
                     file.rstrip(),
                     "--scale-units arcsecperpix --scale-low "+str(pixlo)+" --scale-high "+str(pixhi),
                     "--radius "+str(field),
                     "--downsample "+str(downsample),
                     "--no-plots",#" --no-fits2fits ",
                     "--skip-solved",
                     "--objs "+str(numstars),
                     "--odds-to-tune-up 1e4",
                     "--no-tweak",
                     "--dir",outdir,"--overwrite"])

            rmcmd = "rm -rf "+outdir
            os.system(rmcmd)
            mkdircmd = 'mkdir '+outdir
            os.system(mkdircmd)

            os.system(command)

            outname = fname.split('.')[0]+'.new'
            mvcmd = "mv "+outdir+"/"+outname+" "+datadir+ffinal
            os.system(mvcmd)

            rmcmd = "rm -rf "+outdir
            os.system(rmcmd)

    return

def radec_to_xy(ra,dec,header):
    """
    Overview:
    ---------
    Takes a ra, a dec, and a header and returns the x,y position of the source in the image
    ra and dec are degrees
    x and y are pixel values
    check header for astrometical slove
    """
    if header['CRVAL1']:
        w = wcs.WCS(header)
        world = np.array([[ra, dec]])
        pix = w.wcs_world2pix(world,1) # Pixel coordinates of (RA, DEC)
        x = pix[0,0]
        y = pix[0,1]
        return x,y
    else:
        print(header['OBJECT'] + " has inadequate astrometry information")
        return None,None

def xy_to_radec(x,y,header):
    """
    Overview:
    ---------
    Takes an x, a y, and a header returns ra and dec of source
    x and y are pixel values
    ra and de are degrees
    """
    if header['CRVAL1']:
        w = wcs.WCS(header)
        pix = np.array([[x,y]])
        world = w.wcs_pix2world(pix,1)
        ra = world[0,0]
        dec = world[0,1]
        return ra,dec
    else:
        print(header['OBJECT'] + " has inadequate astrometry information")
        return None,None

#---------------------------------------------------------------------#
# done_in:                                                             #
#----------------------------------------------------------------------#

def done_in(tmaster):

    """
    Overview:
    ---------
    Simple routine to print out the time elapsed since input time

    Calling sequence:
    -----------------
    import time
    tstart = time.time()
    (stuff happens here)
    done_in(tstart)

    """

    t = time.time()
    hour = (t - tmaster)/3600.
    if np.floor(hour) == 1:
        hunit = "hour"
    else:
        hunit = "hours"

    minute = (hour - np.floor(hour))*60.
    if np.floor(minute) == 1:
        munit = "minute"
    else:
        munit = "minutes"

    sec = (minute - np.floor(minute))*60.

    if np.floor(hour) == 0 and np.floor(minute) == 0:
        tout = "done in {0:.2f} seconds"
        print(tout.format(sec))
    elif np.floor(hour) == 0:
        tout = "done in {0:.0f} "+munit+" {1:.2f} seconds"
        print(tout.format(np.floor(minute),sec))
    else:
        tout = "done in {0:.0f} "+hunit+" {1:.0f} "+munit+" {2:.2f} seconds"
        print(tout.format(np.floor(hour),np.floor(minute),sec))

    print(" ")

    return


#----------------------------------------------------------------------#
# get_files:                                                           #
#----------------------------------------------------------------------#

def get_files(d="./",prefix='',tag='',suffix='.fits',clean=True):

    """
    Overview:
    ---------
    Returns list of files with a user defined prefix and suffix withing a
    specified directory


    Calling sequence:
    -----------------
    files = get_files('HATp33b',d='/home/users/bob/stuff/')

    """

    files = glob.glob(d+prefix+"*"+tag+"*"+suffix)

    fct = len(files)

    # Work around for difficult single quote and inconsistent file naming convention
    # due to filter names
    if clean:
        for file in files:
            inname  = file.replace("'","\\'")
            outname =  file.replace("'","")
            if inname != outname:
                mvcmd = "mv "+inname+" "+outname
                os.system(mvcmd)

        files = [file.replace("'","") for file in files]

        for file in files:
            inname  = file
            outname =  file.replace("p.fts",".fts")
            if inname != outname:
                mvcmd = "mv "+inname+" "+outname
                os.system(mvcmd)

        files = [file.replace("p.fts",".fts") for file in files]

    return files,fct




#----------------------------------------------------------------------#
# check_ast                                                            #
#----------------------------------------------------------------------#
def check_ast(file):
    """

    Overview:
    ---------
    Takes an input file and tests to see if there is astrometry in the
    header.

    Calling sequence:
    -----------------
    status = check_ast(file)

    status is 1 if there is no astrometry, else status = 0

    """

    image, header = fits.getdata(file, 0, header=True)
    status = True
    try:
        crval1 = header["CRVAL1"]
    except:
        print("Image has inadequate astrometry information")
        status = False

    return status


#----------------------------------------------------------------------#
# master_bias:
#----------------------------------------------------------------------#

def master_bias(files,write=True,outdir='./',readnoise=False,clobber=False,verbose=True,
                float32=True,tag='',median=False):

    """
   Overview:
    ---------
    Create master bias frame from series of biases (median filter).
    Returns a master_bias frame and writes FITS file to disk in specified
    directory.

    Optionally, the read noise is calculated from the variance of each
    pixel in the bias stack. This is *very* slow. So only use this option
    if you really need to. The readnoise image is also written to disk.

    Inputs:
    -------
    files       : List of flat field files from which a master bias will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)
    readnoise   : Do readnoise calculation (very slow! default False)
    verbose     : Print out progress (default True)

    Calling sequence:
    -----------------
    master_bias = master_bias(biasfiles,write=True,readnoise=False,
                              outdir='/home/users/bob/stuff/')


    """

# Don't redo master_bias unless clobber keyword set
    name  = outdir+'master_bias'+tag+'.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master bias already exists!")
        master_bias = fits.getdata(name,0,header=False)
        return master_bias

# Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz,xsz = image.shape
    stack = np.zeros((fct,ysz,xsz))
    temps = []

# Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = 'Reading {}: frame {} of {} \r'.format(files[i].split('/')[-1],\
                                                          str(i+1),str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        temps.append(header["CCD-TEMP"])
        stack[i,:,:] = image

# Calculate read noise directly from bias frames if prompted
    if readnoise:
        rn = np.zeros((ysz,xsz))
        print("Starting readnoise calculation")
        pbar = tqdm(desc = 'Calculating readnoise', total = ysz, unit = 'rows')
        for i in np.arange(ysz):
            for j in np.arange(xsz):
                rn[i,j] = rb.std(stack[:,i,j])
            pbar.update(1)


# Make a nice plot (after all that hard work)
        aspect = np.float(xsz)/np.float(ysz)
        plt.figure(39,figsize=(5*aspect*1.2,5))
        plt.clf()
        sig = rb.std(rn)
        med = np.median(rn)
        mean = np.mean(rn)
        vmin = med - 2*sig
        vmax = med + 2*sig
        plt.imshow(rn,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
        plt.colorbar()
        plt.annotate(r'$\bar{\sigma}$ = %.2f cts' % mean, [0.95,0.87],horizontalalignment='right',
                     xycoords='axes fraction',fontsize='large')
#                    path_effects=[PathEffects.SimpleLineShadow(linewidth=3,foreground="w")])
        plt.annotate(r'med($\sigma$) = %.2f cts' % med, [0.95,0.8],horizontalalignment='right',
                     xycoords='axes fraction',fontsize='large')
#                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.annotate(r'$\sigma_\sigma$ = %.2f cts' % sig,
                     [0.95,0.73],horizontalalignment='right',
                     xycoords='axes fraction',fontsize='large')
#                    path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
        plt.title("Read Noise")
        plt.xlabel("pixel number")
        plt.ylabel("pixel number")

        if write:
            plt.savefig(outdir+'readnoise'+tag+'.png',dpi=300)

# Calculate master bias frame by median filter
    print('Calculating median of stacked frames...')
    if median:
        master_bias = np.median(stack,axis=0)
    else:
        master_bias = np.mean(stack,axis=0)
        
    # Make a plot
    aspect = np.float(xsz)/np.float(ysz)
    plt.figure(38,figsize=(5*aspect*1.2,5))
    plt.clf()
    sig = rb.std(master_bias)
    med = np.median(master_bias)
    vmin = med - 2*sig
    vmax = med + 2*sig
    plt.imshow(master_bias,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.annotate('Bias Level = %.2f cts' % med, [0.95,0.87],horizontalalignment='right',
                 xycoords='axes fraction',fontsize='large',color='k')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts' % sig, [0.95,0.8],horizontalalignment='right',
                 xycoords='axes fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.95,0.73],horizontalalignment='right',
                 xycoords='axes fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Bias")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

# Write out bias, readnoise and plot
    if write:
        name  = outdir+'master_bias'+tag
        plt.savefig(name+'.png',dpi=300)

        hout = fits.Header()
        hout['CCDTEMP'] = (np.median(temps), "Median CCD temperature")
        hout["TEMPSIG"] = (np.std(temps), "CCD temperature RMS")
        hout["BIAS"] = (med, "Median bias level (cts)")
        hout["BIASSIG"] = (sig, "Bias RMS (cts)")
        if len(glob.glob(name+'.fits')) == 1:
            os.system('rm '+name+'.fits')
        if float32:
            fits.writeto(name+'.fits', np.float32(master_bias), hout)
        else:
            fits.writeto(name+'.fits', master_bias, hout)

        if readnoise:
            name  = outdir+'readnoise'+tag
            if len(glob.glob(name+'.fits')) == 1:
                os.system('rm '+name+'.fits')
            if float32:
                fits.writeto(name+'.fits', np.float32(rn), hout)
            else:
                fits.writeto(name+'.fits', rn, hout)

    return master_bias



#----------------------------------------------------------------------#
# master_dark:
#----------------------------------------------------------------------#

def master_dark(files,bias=None,write=True,outdir='./',clobber=False,float32=True,tag='',
                median=False):
    """
    Overview:
    ---------
    Create master dark frame from series of darks (median filter).
    Returns a master dark frame. If write is specified, a FITS file
    will be written to "outdir" (default is pwd).

    Inputs:
    -------
    files       : List of flat field files from which a master dark will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    bias        : Master bias frame (default None)
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)

    Calling sequence:
    -----------------
    master_dark = master_dark(darkfiles,bias=master_bias,write=True,
                              outdir='/home/users/bob/stuff/')

    """

# Don't redo master_dark unless clobber keyword set
    name  = outdir+'master_dark'+tag+'.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master dark already exists!")
        master_dark = fits.getdata(name,0,header=False)
        return master_dark

 # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    ysz,xsz = image.shape
    stack = np.zeros((fct,ysz,xsz))
    temps = []
    exps = []

# Load stack array and get CCD temperatures
    for i in np.arange(fct):
        output = 'Reading {}: frame {} of {} \r'.format(files[i].split('/')[-1],\
                                                          str(i+1),str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        exp = header["EXPOSURE"]
        exps.append(exp)
        temps.append(header["CCD-TEMP"])
        if length(bias) == 1:
            image = np.float(image)/exp
        else:
            image = (image-bias)/exp
        stack[i,:,:] = image

# Obtain statistics for the master dark image header
    # Temperature
    tmax = np.max(temps)
    tmin = np.min(temps)
    tmean = np.mean(temps)
    tmed = np.median(temps)
    tsig = np.std(temps)
    # Exposure times
    expmax = np.max(exps)
    expmin = np.min(exps)
    print('')
    print("Minimum CCD Temp. %.2f C" % tmin)
    print("Maximum CCD Temp. %.2f C" % tmax)
    print("CCD Temp. rms: %.3f C" % tsig)
    print("CCD Temp. mean: %.2f C" % tmean)
    print("CCD Temp. median: %.2f C" % tmed)

# Create master dark by median filter or mean
    if median:
        master_dark = np.median(stack,axis=0)
    else:
        master_dark = np.mean(stack,axis=0)
        
# Make a plot
    sig = rb.std(master_dark)
    med = np.median(master_dark)
    vmin = med - 2*sig
    vmax = med + 2*sig
    aspect = np.float(xsz)/np.float(ysz)
    plt.figure(37,figsize=(5*aspect*1.2,5))
    plt.clf()
    plt.imshow(master_dark,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.annotate('Dark Current = %.2f cts/sec' % med, [0.72,0.8],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\sigma$ = %.2f cts/sec' % sig, [0.72,0.75],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.annotate(r'$\langle T_{\rm CCD} \rangle$ = %.2f C' % np.median(temps),
                 [0.72,0.7],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
#                 path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])
    plt.title("Master Dark")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")

# Write out plot and master dark array
    if write:
        name = outdir+'master_dark'+tag

        plt.savefig(name+'.png',dpi=300)

        hout = fits.Header()
        hout["TEMPMAX"] = (tmax, "Maximum CCD temperature")
        hout["TEMPMIN"] = (tmin, "Minimum CCD temperature")
        hout["TEMPMED"] = (tmed, "Median CCD temperature")
        hout["TEMPMN"] = (tmean, "Mean CCD temperature")
        hout["TEMPSIG"] = (tsig, "CCD temperature RMS")
        hout["EXPMAX"] = (expmax,"Maximum exposure time")
        hout["EXPMIN"] = (expmin, "Minimum exposure time")
        hout["DARKCNT"] = (med, "Median dark current (cts/sec)")
        hout["DARKSIG"] = (sig, "Dark current RMS (cts/sec)")
        if len(glob.glob(name)) == 1:
            os.system('rm '+name+'.fits')
        if float32:
            fits.writeto(name+'.fits', np.float32(master_dark), hout)
        else:
            fits.writeto(name+'.fits', master_dark, hout)

    return master_dark




#----------------------------------------------------------------------#
# master_flat:
#----------------------------------------------------------------------#

def master_flat(files,bias=None,dark=None,write=True,outdir='./',
                tag='',clobber=False,stretch=3,float32=True,median=False):

    """
    Overview:
    ---------
    Create a master flat using (optionally) a provided bias and dark frame. Output
    is written to "outdir" in FITS format.

    Inputs:
    -------
    files       : List of flat field files from which a master flat will be created.
                  Must be provided, no default.

    Keyword arguments:
    ------------------
    bias        : Master bias frame (default None)
    dark        : Master dark frame calibrated in ADU/sec (default None)
    write       : Toggle to write files to disk (default True)
    outdir      : Directory to which output files are written (default pwd)
    clobber     : Toggle to overwrite files if they already exist in outdir
                  (default False)
    stretch     : Multiple of the noise RMS to stretch image (default 3)


    Calling sequence:
    -----------------
    master_flat = master_flat(flatfiles,bias=master_bias,dark=master_dark,write=True,
                              outdir='/home/users/bob/stuff/')

    """

# Don't redo master_dark unless clobber keyword set
    name = outdir+'master_flat'+tag+'.fits'
    if len(glob.glob(name)) == 1 and not clobber:
        print("Master flat already exists!")
        master_flat = fits.getdata(name,0, header=False)
        return master_flat

 # Get information from inputs and create stack array
    fct = len(files)
    image, header = fits.getdata(files[0], 0, header=True)
    filter = header["filter"]

    ysz,xsz = image.shape
    stack = np.zeros((fct,ysz,xsz))

# Load stack array and get CCD temperatures
    meds = []
    for i in np.arange(fct):
        output = 'Reading {}: frame {} of {} \r'.format(files[i].split('/')[-1],\
                                                          str(i+1),str(fct))
        sys.stdout.write(output)
        sys.stdout.flush()
        image, header = fits.getdata(files[i], 0, header=True)
        image = np.float32(image)
        if header["filter"] != filter:
            sys.exit("Filters do not match!")
        if length(bias) > 1:
            image -= bias
        if length(dark) > 1:
            exptime = header['EXPTIME']
            image -= dark*exptime
        meds.append(np.median(image))
        stack[i,:,:] = image/np.median(image)

# Obtain statistics for the master dark image header
    med = np.median(meds)
    sig = np.std(meds)

# Create master flat by median filter
    if median:
        master_flat = np.median(stack,axis=0)
    else:
        master_flat = np.mean(stack,axis=0)
        
# Make a plot
    sig = rb.std(master_flat)
    med = np.median(master_flat)
    vmin = med - stretch*sig
    vmax = med + stretch*sig
    aspect = np.float(xsz)/np.float(ysz)
    plt.figure(40,figsize=(5*aspect*1.2,5))
    plt.clf()
    plt.imshow(master_flat,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.colorbar()
    plt.title("Master Flat")
    plt.xlabel("pixel number")
    plt.ylabel("pixel number")


# Write out plot and master flat array
    if write:
        plt.savefig(outdir+'master_flat'+tag+'.png',dpi=300)
        hout = fits.Header()
        hout["FILTER"] = (filter, "Filter used when taking image")
        hout["MEDCTS"] = (med, "Median counts in individual flat frames")
        hout["MEDSIG"] = (sig, "Median count RMS in individual flat frames")
        if length(bias) > 1:
            hout.add_comment("Bias subtracted")
        if length(dark) > 1:
            hout.add_comment("Dark subtracted")

        if len(glob.glob(outdir+'master_flat'+tag+'.fits')) == 1:
            os.system('rm '+outdir+'master_flat'+tag+'.fits')
        if float32:
            fits.writeto(outdir+'master_flat'+tag+'.fits', np.float32(master_flat), hout)
        else:
            fits.writeto(outdir+'master_flat'+tag+'.fits', master_flat, hout)

    return master_flat




#----------------------------------------------------------------------#
# brightest_star                                                       #
#----------------------------------------------------------------------#

def brightest_star(image,header,min_distance=10,show=True,threshold_abs=None):

    """
    Under development

    """

    from scipy import ndimage
    from skimage.feature import peak_local_max
    from skimage import data, img_as_float

#    image, header = fits.getdata(file, 0, header=True)
    image = np.float32(image)
    if threshold_abs == None:
        threshold_abs = np.median(image) + 10.0*rb.std(image)

    coordinates = peak_local_max(image, min_distance=min_distance,threshold_abs=threshold_abs)


    peakvals = image[[p[0] for p in coordinates], [p[1] for p in coordinates]]
    peakind = np.argmax(peakvals)
    xpeak = coordinates[peakind,1]
    ypeak = coordinates[peakind,0]


    if show == True:
        plt.figure(1)
        plt.clf()
        sig = rb.std(image)
        med = np.median(image)
        vmin = med - 5*sig
        vmax = med + 15*sig
        plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
        plt.scatter([p[1] for p in coordinates], [p[0] for p in coordinates],marker='o',
                    s=100,facecolor='none',edgecolor='green',linewidth=1.5)

        plt.scatter(xpeak,ypeak,marker='+',s=100,color='green',linewidth=1.5)
        plt.title('Brightest star')

    return np.float64(xpeak),np.float64(ypeak)



#----------------------------------------------------------------------#
# choose_refs:                                                         #
#    routine to select the reference stars to be used
#----------------------------------------------------------------------#

def choose_refs(file,target_ra,target_dec,bias=None,dark=None,flat=None,origin='lower',
                figsize=8,outdir='./',outfile='coordinates.txt',clobber=False):

    '''
    target_ra in hours
    target_dec in degrees
    '''
    
    # Don't redo choose_refs unless clobber keyword set
    if len(glob.glob(outdir+outfile)) == 1 and not clobber:
        print("Reference position file: "+outfile+" already exists!")
        print("... reading saved positions")
        coords = np.loadtxt(outdir+outfile)
        ras = coords[:,0]
        decs = coords[:,1]
        return ras,decs

    # Convert ra and dec strings to decimal degree floats
    coords = SkyCoord(target_ra,target_dec,unit=(u.hour,u.deg))
    ra0  = coords.ra.deg
    dec0 = coords.dec.deg

    # Read image
    image,header = fits.getdata(file, 0, header=True)
    image = image.astype('float')
    ysz,xsz = image.shape
    if length(bias) > 1:
        image -= bias
    if length(dark) > 1:
        image -= (dark*header['EXPTIME'])
    if length(flat) > 1:
        image /= flat

# Read header
    hdulist = astropy.io.fits.open(file)

# Convert RA and Dec to pixel coordinates
    w = wcs.WCS(hdulist[0].header)
    world0 = np.array([[ra0, dec0]])
    pix0 = w.wcs_world2pix(world0,1) # Pixel coordinates of (RA, DEC)
    x0 = pix0[0,0]
    y0 = pix0[0,1]

    if (x0 > xsz) |(x0 < 0) | (y0 > ysz) | (y0 < 0):
        print x0,y0
        print("Target star position is out of range!")
        return None, None

# Get image information
    sig = rb.std(image)
    med = np.median(image)
    vmin = med - 5*sig
    vmax = med + 15*sig

# Plot image
    plt.ion()
    fig = plt.figure(99,figsize=(figsize,figsize))
    ax = fig.add_subplot(111)
    ax.imshow(image,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin=origin)
    ax.scatter(x0,y0,marker='o',s=100,facecolor='none',edgecolor='green',linewidth=1.5)
    ax.scatter(x0,y0,marker='+',s=100,facecolor='none',edgecolor='green',linewidth=1.5)
    plt.draw()

# Click image and receive x and y values
    refx=[]
    refy=[]
    def onclick(event):
        newx = event.xdata
        newy = event.ydata
        if newx < xsz and newy < ysz and newx > 0 and newy > 0:
            refx.append(newx)
            refy.append(newy)
            print("------------------------------")
            print("refx = " , newx)
            print("refy = " , newy)
            ax.scatter(refx,refy,marker='o',s=100,facecolor='none',edgecolor='yellow',linewidth=1.5)
            plt.draw()

# Engage "onclick"
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

# Stall here such that canvas can disconnect before further calculations
    print("Click positions of reference stars")
    print(raw_input("Press return when finished selecting sources \n"))

# Disengage "onclick"
    fig.canvas.mpl_disconnect(cid)

# Convert xs and ys to RAs and Decs
    raval,decval = w.wcs_pix2world(refx,refy,1)


# Put the target star at the beginning of the lists
    ras = np.append(ra0,raval)
    decs = np.append(dec0,decval)

# Write out coordinate file
    coords = np.zeros((len(ras),2))
    coords[:,0] = ras
    coords[:,1] = decs
    np.savetxt(outdir+outfile,coords)

    return ras,decs




#----------------------------------------------------------------------#
# choose_target:                                                       #
#    routine to select a target from a file                            #
#----------------------------------------------------------------------#

def choose_target(file,target_ra,target_dec,bias=None,dark=None,flat=None,
                   figsize=8):

# Convert ra and dec strings to decimal degree floats
    coords = SkyCoord(target_ra,target_dec,unit=(u.deg,u.deg))
    ra0  = coords.ra.deg
    dec0 = coords.dec.deg

# Read image
    image = fits.getdata(file, 0, header=False).astype('float')
    ysz,xsz = image.shape

    if length(bias) > 1:
        image -= bias
    if length(dark) > 1:
        image -= dark
    if length(flat) > 1:
        image /= flat

# Read header
    hdulist = fits.open(file)

# Convert RA and Dec to pixel coordinates
    w = wcs.WCS(hdulist[0].header)
    world0 = np.array([[ra0, dec0]])
    pix0 = w.wcs_world2pix(world0,1) # Pixel coordinates of (RA, DEC)
    x0 = pix0[0,0]
    y0 = pix0[0,1]

    if (x0 > xsz) |(x0 < 0) | (y0 > ysz) | (y0 < 0):
        print("Target star position is out of range!")
        return

    # Get image information
    sig = rb.std(image)
    med = np.median(image)
    vmin = med - 5*sig
    vmax = med + 15*sig

    # Click image and receive x and y values
    clickvals = {}
    def onclick(event):
        xclick = event.xdata
        yclick = event.ydata
        clickvals['x'] = xclick
        clickvals['y'] = yclick
        #xref = np.append(xref,xclick)
        #yref = np.append(yref,yclick)
        if xclick < xsz and yclick < ysz and xclick > 0 and yclick > 0:
            print("------------------------------")
            print("x = " , xclick)
            print("y = " , yclick)
            ax.scatter(xclick,yclick,marker='o',s=100,facecolor='none',edgecolor='yellow',linewidth=1.5)
            plt.draw()

    # Plot image
    plt.ioff()
    fig = plt.figure(99,figsize=(figsize,figsize))
    ax = fig.add_subplot(111)
    ax.imshow(image,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    ax.scatter(x0,y0,marker='o',s=100,facecolor='none',edgecolor='green',linewidth=1.5)
    ax.scatter(x0,y0,marker='+',s=100,facecolor='none',edgecolor='green',linewidth=1.5)
    plt.draw()

    # Engage "onclick"
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.waitforbuttonpress()

    xclick = clickvals['x']
    yclick = clickvals['y']

    # Disengage "onclick"
    fig.canvas.mpl_disconnect(cid)

    # Convert xs and ys to RAs and Decs
    raval,decval = w.wcs_pix2world(xclick,yclick,1)

    return raval,decval




def check_skyrad(file=None,image=None,header=None,ra=None,dec=None,bias=None,dark=None,flat=None,
                 skyrad=[20,25],aperture=None,figsize=(8,8),siglo=3.0,sighi=5.0):

    if file is not None:
        status = check_ast(file)
        if status == 0:
            print 'No astrometry in header'
            return

        image, header = fits.getdata(file, 0, header=True)
        image = np.float32(image)
    elif image is not None:
        image = np.float32(image)
    else:
        print 'File or image must be supplied'
        return
    
    if np.shape(bias) == np.shape(image):
        image -= bias
    if np.shape(dark) == np.shape(image):
        image -= (dark*header['EXPTIME'])
    if np.shape(flat) == np.shape(image):
        image /= flat

    dict = optimal_aperture(image,header,ra=ra,dec=dec,skyrad=skyrad,plot=False)
    x = dict['xcen']
    y = dict['ycen']
    if not aperture:
        aperture = dict['optimal_aperture']

    sz = int(round(max(200,np.max(skyrad)*2.25)))

    # Check this code !!!
    if sz % 2 == 0:
        sz += 1

    yround = int(round(y))
    xround = int(round(x))
    patch = image[yround-sz/2:yround+sz/2+1,xround-sz/2:xround+sz/2+1]

    position = (x-np.round(x)+sz/2,y-np.round(y)+sz/2)
    
    bkg_aperture = CircularAnnulus(position, r_in=skyrad[0], r_out=skyrad[1])
    aperture = CircularAperture(position, r=aperture)

    sig = rb.std(patch)
    med = np.median(patch)
    vmin = med - siglo*sig
    vmax = med + sighi*sig
    plt.figure(20,figsize=figsize)
    plt.clf()
    plt.imshow(patch,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.scatter(sz/2,sz/2,marker='+',s=200,color='yellow',linewidth=1.5,label='Original position')
    plt.title('zoom in of target')
    aperture.plot(color='cyan')
    bkg_aperture.plot(color='cyan', hatch='//', alpha=0.8)

    # Need to see if star is saturated!!
    # This may be broken
    #xpos = np.arange(-(sz/2),(sz/2)+1,1)
    #ypos = np.arange(-(sz/2),(sz/2)+1,1)
    xpos = np.arange(-(20),(20)+1,1)
    ypos = np.arange(-(20),(20)+1,1)
    X,Y = np.meshgrid(xpos,ypos)
    fig = plt.figure(21,figsize=figsize)
    ax = fig.gca(projection='3d')
    ax.set_title('Point Spread Function')
    zoom = patch[sz/2-20:sz/2+21,sz/2-20:sz/2+21]
    surf = ax.plot_surface(X, Y, zoom, rcount=100, ccount=100,
                           cmap=cm.gist_earth_r)
#    ax.set_xlim(-20,20)
#    ax.set_ylim(-20,20)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    return



#----------------------------------------------------------------------#
# optimal_aperture:                                                    #
#            calculate optimal aperture for a source at "x" and "y"    #
#            coords in "image" using defined skyradii = [in,out]       #
#----------------------------------------------------------------------#

def optimal_aperture(image,header,x=None,y=None,ra=None,dec=None,
                     bias=None, dark=None, flat=None, apmax=None,
                     aperture=None,skyrad=[20,25],use_old=False,
                     plot=False):

    """
    optimal_aperture takes an image, a header, and skyrad [in,out]
        either x and y in pixel values, or ra and dec in degrees
    returns dictionary at end of function
    """

    if ra == None and dec == None:
        ra,dec = xy_to_radec(x,y,header)
    elif x == None and y == None:
        x,y = radec_to_xy(ra,dec,header)

    exptime = header['EXPTIME']
    if length(bias) > 1:
        image -= bias
    if length(dark) > 1:
        image -= (dark*exptime)
    if length(flat) > 1:
        image /= flat
        
    jd = header["jd"] + (header["exptime"]/2.0)/(24.0*3600.0)
    
    try:
        secz = header['airmass']

    except:
        obs.date = ephem.date(jd-2415020.0) # pyephem uses Dublin Julian Day (why!!??!?)
        star = ephem.FixedBody()

        star._ra = np.radians(ra)
        star._dec = np.radians(dec)
        star.compute(obs)
        secz = airmass(np.degrees(star.alt))

    # Create zoom in of target
    sz = 60
    try:
        patch = image[int(round(y-sz/2)):int(round(y+sz/2)),int(round(x-sz/2)):int(round(x+sz/2))]
    except:
        dict = {'optimal_aperture':[], 'optimal_flux':[],
            'optimal_fluxerr':[], 'xcen':[], 'ycen':[],
            'fwhm':[], 'aspect':[],'snrmax':[],'totflux':[],
            'totfluxerr':[], 'totflux_aperture':[],'chisq':[],
                'curve_of_growth':[],'secz':[], 'jd':[]}
        return dict

# Plot zoom in with fits overlaid
    sig = rb.std(patch)
    med = np.median(patch)
    vmin = med - 5*sig
    vmax = med + 15*sig
    if plot:
        plt.figure(2)
        plt.clf()
        plt.imshow(patch,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
        plt.scatter(sz/2,sz/2,marker='+',s=200,color='yellow',linewidth=1.5,label='Original position')
        plt.title('zoom in of target star')

    xg = np.linspace(0, sz-1, sz)
    yg = np.linspace(0, sz-1, sz)
    xg, yg = np.meshgrid(xg, yg)

    #mom = moments(patch)
    #p0 = [mom[0],mom[1],mom[2],mom[3],mom[4],0,mom[5]]

    # For starting parameters assume seeing ~ 4"
    # p = [height,x,y,width_x,width_y,angle,baselevel]
    sigguess = 4/2.355/0.61
    p0 = [np.max(patch[int(round(sz/2-sigguess)):int(round(sz/2+sigguess)),
                       int(round(sz/2-sigguess)):int(round(sz/2+sigguess))]),
          sz/2,sz/2,sigguess,sigguess,0,med]
    

# Fit 2D Guassian to target
    try:
        params, pcov = opt.curve_fit(Gaussian2D, (xg, yg), patch.ravel(), p0=p0)
        base = params[-1]
        peak = params[0]
        fwhm = np.sqrt(np.abs(params[3]*params[4]))*2.0*np.sqrt(2*np.log(2))
        aspect = min(np.abs(params[3]),np.abs(params[4]))/max(np.abs(params[3]),np.abs(params[4]))
#        norm = params[0] * 2.0 * np.pi * np.abs(params[3] * params[4])
        level = peak*np.array([0.1, 1./np.e, 0.95]) + params[-1]
        xpeak = params[1] + x - sz/2
        ypeak = params[2] + y - sz/2

        if plot:
            fit = Gaussian2D((xg, yg), *params)
            plt.contour(xg,yg,fit.reshape(sz,sz),level,colors='cyan')
            plt.scatter(params[1],params[2],marker='+',s=200,color='cyan',linewidth=1.5,label='Updated position')
    except:
        print('Fit Failed')
        aspect = np.nan ; fwhm = np.nan ; norm = np.nan ; peak = np.nan
        fit = np.nan ; level = np.nan
        xpeak = x ; ypeak = y


# Create vector of apertures
    if not apmax:
        apmax = np.min(skyrad)-1
    ap = np.linspace(1,apmax,100)


    '''
    #Using photutils
    position = [(xpeak,ypeak)]
    bkg_apertures = CircularAnnulus(position, r_in=skyrad[0], r_out=skyrad[1])
    flux = []
    
    for val in ap:
        circ_ap = CircularAperture((xpeak,ypeak), r=val)
        phot = aperture_photometry(image, circ_ap,  method='exact')
        bkg = aperture_photometry(image, bkg_apertures)
        bkg_mean = bkg['aperture_sum'] / bkg_apertures.area()
        bkg_sum = bkg_mean * circ_ap.area()
        flux_bkgsub = phot['aperture_sum'] - bkg_sum
        phot['bkg_sum'] = bkg_sum
        phot['aperture_sum_bkgsub'] = flux_bkgsub
    xval = xpeak
    yval = ypeak
    dy = 0
    dx = 0

    pdb.set_trace()

    '''


# y and x are swapped due to python vs. IDL indexing !!!
# (could also pass image.T)
#    phot = djs.djs_phot(y,x,ap,skyrad,image,skyrms=True,flerr=True,skyval=True,cbox=10)
    phot = djs.djs_phot(ypeak,xpeak,ap,skyrad,image,skyrms=True,flerr=True,skyval=True,cbox=10)
    xval = phot["xcen"]
    yval = phot["ycen"]
    dy = yval-xpeak
    dx = xval-ypeak
    
# Chi Squared of Gaussian fit
    patchrms = np.sqrt(patch)
    xi,yi = np.where(patchrms == 0.0)
    if len(xi) > 0:
        patchrms[xi,yi] = 1.0
    try:
        chisq = np.sum((patch-fit.reshape(sz,sz))**2/patchrms**2)/(sz**2 - len(params) - 1.0)
    except:
        chisq =  np.nan

# Optimize based on signal to noise from counts and background RMS
    flux = phot["flux"][:,0]
    skyrms = phot["skyrms"][0]
    fluxerr = phot['fluxerr'][:,0]
    sbg = skyrms*ap*np.sqrt(np.pi) # rms * sqrt(# of pixels in phot aperture)
    snr = flux/np.sqrt(sbg**2 + flux)
    snrmax = np.max(snr)
    maxi = np.argmax(flux)
    totflux = flux[maxi]
    totfluxerr = fluxerr[maxi]
    totap = ap[np.argmax(flux)]

# Plot optimization
    cog = flux/totflux
    if plot:
        plt.figure(3)
        plt.clf()
        plt.subplot(2,1,1)
        plt.plot(ap,cog)
        plt.xlabel("aperture radius (pixels)")
        plt.ylabel("normalized flux")
        plt.annotate('Total Counts = %.f' % np.max(flux), [0.86,0.57],horizontalalignment='right',
                     xycoords='figure fraction',fontsize='large')
        plt.ylim([0,1.1])
        plt.axhline(y=1.0,linestyle='--',color='green')
        plt.subplot(2,1,2)
        plt.plot(ap,snr)
        plt.xlabel("aperture radius (pixels)")
        plt.ylabel("SNR")

# Optimal aperture
    op_ap = ap[np.argmax(snr)]
    op_ap_flux = flux[np.argmax(snr)]
    op_ap_fluxerr = fluxerr[np.argmax(snr)]

    if plot:
        plt.axvline(x=op_ap,linestyle='--',color='red')
        plt.annotate('SNR maximum = %.f' % np.max(snr), [0.8,0.17],horizontalalignment='right',
                     xycoords='figure fraction',fontsize='large')
        plt.annotate('Optimal aperture = %.2f' % op_ap, [0.8,0.12],horizontalalignment='right',
                     xycoords='figure fraction',fontsize='large')
        if aperture != None:
            plt.axvline(x=aperture,linestyle='--',color='black')

        plt.draw()

    dict = {'optimal_aperture':op_ap, 'optimal_flux':op_ap_flux,
            'optimal_fluxerr':op_ap_fluxerr, 'xcen':yval, 'ycen':xval,
            'fwhm':fwhm, 'aspect':aspect,'snrmax':snrmax,'totflux':totflux,
            'totfluxerr':totfluxerr, 'totflux_aperture':totap,'chisq':chisq,
            'curve_of_growth':[ap,cog],'secz':secz,'jd':jd, 'exptime':exptime}

    if use_old:
        return op_ap,xval,yval,fwhm,aspect,snrmax,totflux,totap,chisq,secz
    else:
        return dict



#----------------------------------------------------------------------#
# total_flux
#----------------------------------------------------------------------#
def total_flux(file,obsc=0.47,gain=1.7,SpT=None,skyrad=[30,40],mag=None,
               dark=None,bias=None,flat=None,ra=None,dec=None,
               lat='34 28 00.5',lon='-119 10 38.5',elevation=494.7,
               doplot=False,outdir='./',name='integrand',
               camera='U16m',filter='V',brightest=False,nearest=True,
               network='swift'):

    import ephem

    dict = { "secz":[], "flux":[], "fluxerr":[], "tau":[], "tauerr":[],
            "exptime":[], "jd":[], "snr":[], "chisq":[], "aperture":[],
             "xpos":[], "ypos":[]}


    dpath, stpath = get_info_path()
    dpath += '/'
    stpath += '/'

    obs = ephem.Observer()
    obs.lat = lat
    obs.lon = lon
    obs.elevation=elevation

    if mag == None:
        sys.exit("Must supply source magnitude")
    if SpT == None:
        sys.exit("Must supply source spectral type")

    image,header = fits.getdata(file, 0, header=True)
    image = np.float32(image)
    exptime = header["exptime"]
    if ra == None:
        object = header["object"]
        object.replace('_',' ')
        targ = SkyCoord.from_name(object)
        ra = targ.ra.degree
        dec = targ.dec.degree

    jd = header["jd"] + (exptime/2.0)/(24.0*3600.0)
    dict["jd"] = jd
    obs.date = ephem.date(jd-2415020.0) # pyephem uses dublin julian day (why!!??!?)

    coord = SkyCoord(ra, dec, unit='deg')
    h = coord.to_string('hmsdms').split('s')
    print h
    loc = [] # location (ra, dec)
    for s in h[:2]:
        s = s.replace(' ', '')
        if s[0]=='-' or s[0]=='+':
            loc.append(s[1:3]+':'+s[4:6]+':'+s[7:])
        else:
            loc.append(s[:2]+':'+s[3:5]+':'+s[6:])

    star = ephem.FixedBody()
    # Not sure this will work!
    star._ra = loc[0]#np.radians(ra)
    star._dec = loc[1]#de
    star.compute(obs)
#    dict["secz"] = 1.0/np.cos(np.pi/2.0 - star.alt)
    dict['secz'] = header['airmass']

    if brightest == True:
        xpeak,ypeak = brightest_star(image,header)

    elif nearest == True:
        hdulist = astropy.io.fits.open(file)
        w = wcs.WCS(hdulist[0].header)
        world0 = np.array([[ra, dec]])
        pix0 = w.wcs_world2pix(world0,1) # Pixel coordinates of (RA, DEC)
        #print 'calc secz:',dict['secz']
        xpeak = pix0[0,0]
        ypeak = pix0[0,1]
        #plt.ion()
        #plt.figure(934)
        sz = 100
        patch = image[int(round(ypeak-sz/2)):int(round(ypeak+sz/2)),int(round(xpeak-sz/2)):int(round(xpeak+sz/2))]
        xg = np.linspace(0, sz-1, sz)
        yg = np.linspace(0, sz-1, sz)
        xg, yg = np.meshgrid(xg, yg)

        sigguess = 4/2.355/0.61
        p0 = [np.max(patch[int(round(sz/2-sigguess)):int(round(sz/2+sigguess)),
                           int(round(sz/2-sigguess)):int(round(sz/2+sigguess))]),
              sz/2,sz/2,sigguess,sigguess,0,med]

        #mom = moments(patch)
        #p0 = [mom[0],mom[1],mom[2],mom[3],mom[4],0,mom[5]]
        # Try alternative fitting!
        params, pcov = opt.curve_fit(Gaussian2D, (xg, yg), patch.ravel(), p0=p0)
        base = params[-1]
        peak = params[0]
        fwhm = np.sqrt(np.abs(params[3]*params[4]))*2.0*np.sqrt(2*np.log(2))
        aspect = min(np.abs(params[3]),np.abs(params[4]))/max(np.abs(params[3]),np.abs(params[4]))
        xpeak = params[2] + xpeak - sz/2
        ypeak = params[1] + ypeak - sz/2

        #plt.imshow(patch,vmin=np.median(patch)-2*np.std(patch),vmax=np.median(patch)+5*np.std(patch))
#        params = fitgaussian(patch)
#        fit = gaussian(*params)
#        level = (np.max(patch) - np.median(patch))*np.array([0.95,0.5,0.1])
        #plt.contour(fit(*indices(patch.shape)),level,colors='blue')
#        aspect = min(params[3],params[4])/max(params[3],params[4])
#        fwhm = np.sqrt(params[3]*params[4])*2.0*np.sqrt(2*np.log(2))
#        norm = params[0] * 2.0 * np.pi * params[3] * params[4]
#        peak = params[0]
#        xpeak = params[2] + xpeak - sz/2
#        ypeak = params[1] + ypeak - sz/2

    else:
        hdulist = astropy.io.fits.open(file)
        w = wcs.WCS(hdulist[0].header)
        world0 = np.array([[ra, dec]])
        pix0 = w.wcs_world2pix(world0,1) # Pixel coordinates of (RA, DEC)
        xpeak = pix0[0,0]
        ypeak = pix0[0,1]


    if length(bias) > 1:
        image -= bias
    if length(dark) > 1:
        image -= (dark*exptime)
    if length(flat) > 1:
        image /= flat

    opdict = optimal_aperture(image, header, x=xpeak,y=ypeak,skyrad=skyrad)

    #ap,xval,yval,fwhm,aspect,snrmax,totflux,totap,chisq

    dict["snr"] = opdict['snrmax']
    dict["fwhm"] = opdict['fwhm']
    dict["aspect"] = opdict['aspect']


    phot = djs.djs_phot(ypeak,xpeak,opdict['totflux_aperture'],skyrad,image,
                        skyrms=True,flerr=True,skyval=True,cbox=10)

    flux = phot["flux"][0]/exptime
    fluxerr = phot["fluxerr"][0]/exptime
    dict["xpos"] = xpeak
    dict["ypos"] = ypeak

    Aeff = np.pi*( (70.0 / 2.0)**2 - (70.0 * obsc /2.0)**2 )


    lam = np.linspace(4250,7000,endpoint=True,num=1000)

    if camera == 'U16m':
        QEdata = np.loadtxt(dpath+'QE_U16m.txt')
        QEfunc = interp(QEdata[:,0]*10.0,QEdata[:,1]/100.0,kind='cubic')
        QEint = QEfunc(lam)
    elif camera == 'STL11000':
        QEdata = np.loadtxt(dpath+'QE_STL11000.txt')
        QEfunc = interp(QEdata[:,0]*10.0,QEdata[:,1]/100.0,kind='cubic')
        QEint = QEfunc(lam)
    elif camera == 'STL6303':
        QEdata = np.loadtxt(dpath+'QE_STL6303.txt')
        QEfunc = interp(QEdata[:,0]*10.0,QEdata[:,1]/100.0,kind='linear')
        QEint = QEfunc(lam)
    elif camera == 'iKON-L':
        QEdata = np.loadtxt(dpath+'QE_iKON-L.txt')
        QEfunc = interp(QEdata[:,0]*10.0,QEdata[:,1]/100.0,kind='cubic')
        QEint = QEfunc(lam)

    frefdata = np.loadtxt(dpath+'Vpassband.dat')
    freffunc = interp(frefdata[:,0],frefdata[:,1],kind='cubic')
    frefint = freffunc(lam)
    inds = np.where(frefint < 0)
    frefint[inds] = 0.0

    if filter == 'V':
        fdata = np.loadtxt(dpath+'AstrodonTransmissionCurves.txt',skiprows=2)
        ffunc = interp(fdata[:,0]*10.0,fdata[:,3]/100,kind='linear')
        fint = ffunc(lam)
    if filter == 'G':
        fdata = np.loadtxt(dpath+'BaaderG.txt')
        ffunc = interp(fdata[:,0]*10.0,fdata[:,1]/100,kind='linear')
        fint = ffunc(lam)
    elif filter == 'Lum':
        fdata = np.loadtxt(dpath+'Lumpassband.txt')
        ffunc = interp(fdata[:,0]*10.0,fdata[:,1]/100,kind='linear')
        fint = ffunc(lam)


    sdata = np.loadtxt(stpath+SpT+'.dat')
    starfunc = interp(sdata[:,0],sdata[:,1],kind='linear')
    starint = starfunc(lam)
    starflux = starfunc(5500.0)

    refdata = np.loadtxt(stpath+'a0v.dat')
    reffunc = interp(refdata[:,0],refdata[:,1],kind='linear')
    refint = reffunc(lam)
    calflux = reffunc(5500.0)


# Make a nice plot
    if doplot:
        plt.figure(11)
        plt.clf()
        plt.plot(lam,fint,'-b',label=r'$'+filter+'$ (used)')
        plt.plot(lam,frefint,'-k',label=r'$V$ (standard)')
        plt.plot(lam,QEint,'-c',label='Q.E.')
        plt.plot(lam,0.5*refint/calflux,color='gray',label='A0V (standard)')
        plt.plot(lam,0.5*starint/starflux,color='purple',label='Stellar spectrum')
        plt.axvline(x=5500,color='purple',linestyle='--')
        arg = fint*0.5*starint/calflux*QEint
        plt.fill_between(lam,arg,hatch='///',label='Integrand',color="none",edgecolor='orange')
        plt.legend(loc=1)
        plt.xlim([np.min(lam),np.max(lam)])
        plt.xlabel(r'Wavelength ($\AA$)')
        plt.ylabel('Transmission')
        plt.savefig(outdir+name+'_int.png',dpi=300)


    mf =  sp.integrate.simps(frefint*refint,x=lam)/sp.integrate.simps(frefint,x=lam)

# Zeropoint flux of A0V star in Johnson V = 3.631e-9 ergs/cm2/s/A
# Bessell et al. 1998
    norm = 3.631e-9 / mf
    starint = starint * norm

    tau  = flux * 10**(0.4 * mag) * c.h *c.c * gain / \
        (Aeff * sp.integrate.simps(fint*starint*QEint*lam/1e8,x=lam))

    tauerr = fluxerr * 10**(0.4 * mag) * c.h *c.c * gain / \
        (Aeff * sp.integrate.simps(fint*starint*QEint*lam/1e8,x=lam))

    dict["flux"] = flux
    dict["fluxerr"] = fluxerr
    dict["tau"] = tau
    dict["tauerr"] = tauerr
    dict["exptime"] = exptime
    dict["chisq"] = opdict['chisq']
    dict["aperture"] = opdict['totflux_aperture']

    return dict



def batch_total_flux(files,ra=None,dec=None,mag=None,SpT=None,flat=None,object=None,dark=None,bias=None,
                     outdir='./',filter='V',camera='U16m',gain=None,brightest=False,nearest=True,
                     skyrad=[30,40],network='swift'):

    if camera == 'U16m' and gain == None:
        # gain = 1.4 (measured?)
        gain = 1.7 # (spec'ed)
    elif camera == 'STL11000' and gain == None:
        gain = 0.77 # (measured?)
    elif camera == 'STL6303' and gain == None:
        gain = 2.3 # (spec'ed)
    elif camera == 'iKON-L' and gain == None:
        gain = 8.1

    print("Using a camera gain of %.2f" % gain)

    if mag == None:
        sys.exit("Must supply source magnitude")
    if SpT == None:
        sys.exit("Must supply source spectral type")

    image,header = fits.getdata(files[0], 0, header=True)
    image = np.float32(image)


    if not object:
        object = header["object"]

    if ra == None or dec == None:
        #obj = object.replace('_',' ')
        targ = SkyCoord.from_name(object)
        ra = targ.ra.degree
        dec = targ.dec.degree


    data = {"tau":[], "tauerr":[], "flux":[], "fluxerr":[], "secz":[],
            "fwhm":[], "aspect":[], "chisq":[], "jd":[], "aperture":[],
            "xpos":[], "ypos":[], "exptime":[]}

    for i in np.arange(len(files)):
        doplot = i == 0
        print("File: "+files[i])
        info = total_flux(files[i],mag=mag,SpT=SpT,ra=ra,dec=dec,name=object,gain=gain, \
                          doplot=doplot,outdir=outdir,filter=filter,camera=camera, \
                          brightest=brightest,flat=flat,bias=bias,dark=dark,skyrad=skyrad,
                          network=network)
        data["tau"]      = np.append(data["tau"],    info["tau"])
        data["tauerr"]   = np.append(data["tauerr"], info["tauerr"])
        data["flux"]     = np.append(data["flux"],   info["flux"])
        data["fluxerr"]  = np.append(data["fluxerr"],info["fluxerr"])
        data["secz"]     = np.append(data["secz"],   info["secz"])
        data["fwhm"]     = np.append(data["fwhm"],   info["fwhm"])
        data["aspect"]   = np.append(data["aspect"], info["aspect"])
        data["chisq"]    = np.append(data["chisq"],  info["chisq"])
        data["jd"]       = np.append(data["jd"],     info["jd"])
        data["aperture"] = np.append(data["aperture"],info["aperture"])
        data["xpos"]     = np.append(data["xpos"],   info["xpos"])
        data["ypos"]     = np.append(data["ypos"],   info["ypos"])
        data["exptime"]  = np.append(data["exptime"],info["exptime"])

    return data


def exponential(params, *args):
    x = args[0]
    y = args[1]
    yerr = args[2]
    amp, index = params
    y_model = amp*np.exp(index*x)
    error = (y-y_model)/yerr
    return np.sum(error**2)


#----------------------------------------------------------------------#
# do_phot:                                                             #
#    given a file and the ra and dec coordinates of stars, return
#    vectors of JD, flux, flux err, sky and skyrms
#----------------------------------------------------------------------#


def do_phot(file=None,ras=None,decs=None,image=None,header=None,
            aperture=None,skyrad=np.array([32,42]),dark=None,peak=None,
            bias=None,flat=None,cr_reject=False,plot=False,figsize=(10,10)):
    '''

    '''
# Get image and header
    if type(file) != types.NoneType:
        image, header = fits.getdata(file, 0, header=True)
        image = np.float32(image)
    elif type(image) ==  types.NoneType or type(header) ==  types.NoneType:
        print("You must supply a filename or an image and header")
        return None
        
    exptime = header['EXPTIME']

    cd = np.array([[header['CD1_1'],header['CD1_2']],
                   [header['CD2_1'],header['CD2_2']]])

    plate_scale = np.mean(np.abs(np.linalg.eig(cd)[0])*3600.0)
    
# Apply calibrations if provided
    if length(bias) > 1:
        image -= bias
    if length(dark) > 1:
        image -= (dark*exptime)
    if length(flat) > 1:
        image /= flat

    if cr_reject:
        image, mask = cosmicray_lacosmic(image, sigclip=5)
        
# Get image info
    jd = header["JD"]
    image = np.array(image)
    ysz, xsz = np.shape(image)

# Get image astrometry
    #hdulist = astropy.io.fits.open(file)
    #w = wcs.WCS(hdulist['PRIMARY'].header)
    w = wcs.WCS(header)
    xpix,ypix = w.wcs_world2pix(ras,decs,1) # Pixel coordinates of (RA, DEC)

# Image characteristics and plot
    if plot:
        sig = rb.std(image)
        med = np.median(image)
        vmin = med - 3*sig
        vmax = med + 5*sig
        plt.figure(1,figsize=figsize)
        plt.clf()
        plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
        plt.scatter(xpix,ypix,marker='o',s=100,facecolor='none',edgecolor='yellow',linewidth=1.5)
        plt.scatter(xpix[0],ypix[0],marker='+',s=100,facecolor='none',edgecolor='green',linewidth=1.5)
        plt.title('target field')
        plt.show()

    # Get the optimal aperture (max SNR) for the image based on the target star
    #dict = optimal_aperture(image,header,x=xpix[0],y=ypix[0],skyrad=skyrad,aperture=aperture)

    if aperture != None:
        ap = aperture
    else:
        ap = dict['optimal_aperture']
        print 'colin and piper go fix your code'

    if length(skyrad) == 2:
        skyrad = np.tile(skyrad,(length(ras),1))

# Loop to catch reference sources that may fall off chip
    flux   = np.zeros(len(xpix))
    flerr  = np.zeros(len(xpix))
    skyval = np.zeros(len(xpix))
    skyrms = np.zeros(len(xpix))
    peaks  = np.zeros(len(xpix))
    for i in np.arange(len(xpix)):
        buf = 10
        xtest = xpix[i]
        ytest = ypix[i]
        sr = skyrad[i]
        if (xtest >= buf) & (xtest < xsz-buf) & (ytest >= buf) & (ytest < ysz-buf):
            phot = djs.djs_phot(ytest,xtest,ap,sr,image,skyrms=True,flerr=True,skyval=True,cbox=10,peak=True)
            flux[i]   = phot["flux"]
            flerr[i]  = phot["fluxerr"]
            skyval[i] = phot["skyval"]
            skyrms[i] = phot["skyrms"]
            peaks[i] = phot["peak"]
        else:
            print "star not in image!"

    # Append to the optimal aperture dictionary and return
    dict['jd'] = jd
    dict['exptime'] = exptime
    dict['flux'] = flux
    dict['flux_err'] = flerr
    dict['skyval'] = skyval
    dict['skyrms'] = skyrms
    dict['aperture'] = ap
    dict['plate_scale'] = plate_scale
    dict['peak'] = peaks
    
    return dict

#----------------------------------------------------------------------#
# simple_phot:                                                             #
#    given a file and the ra and dec coordinates of stars, return
#    vectors of JD, flux, flux err, sky and skyrms
#----------------------------------------------------------------------#


def simple_phot(file,ras,decs,aperture=None,skyrad=np.array([32,42]),dark=None,bias=None,flat=None):

# Get image and header
    image, header = fits.getdata(file, 0, header=True)

# Apply calibrations if provided
    if length(bias) > 1:
        image -= bias
    if length(dark) > 1:
        image -= dark
    if length(flat) > 1:
        image /= flat

# Get image info
    jd = header["JD"]
    exptime = header["EXPTIME"]
    image = np.array(image)
    ysz, xsz = shape(image)

# Get image astrometry
    hdulist = astropy.io.fits.open(file)
    w = wcs.WCS(hdulist['PRIMARY'].header)
    xpix,ypix = w.wcs_world2pix(ras,decs,1) # Pixel coordinates of (RA, DEC)

# Image characteristics and plot
    sig = rb.std(image)
    med = np.median(image)
    vmin = med - 3*sig
    vmax = med + 5*sig
    plt.figure(1)
    plt.clf()
    plt.imshow(image,vmin=vmin,vmax=vmax,cmap='gist_heat',interpolation='nearest',origin='lower')
    plt.scatter(xpix,ypix,marker='o',s=100,facecolor='none',edgecolor='yellow',linewidth=1.5)
    plt.scatter(xpix[0],ypix[0],marker='+',s=100,facecolor='none',edgecolor='green',linewidth=1.5)
    plt.title('target field')
    plt.draw()

# Get the optimal aperture (max SNR) for the image based on the target star
    ap,xval,yval,fwhm,aspect,snr,totflux,totap,chisq = \
        optimal_aperture(xpix[0],ypix[0],image,skyrad,aperture=aperture)

    if aperture != None: ap = aperture

# Loop to catch reference sources that may fall off chip
    flux   = np.zeros(len(xpix))
    flerr  = np.zeros(len(xpix))
    skyval = np.zeros(len(xpix))
    skyrms = np.zeros(len(xpix))
    for i in np.arange(len(xpix)):
        buf = 10
        xtest = xpix[i]
        ytest = ypix[i]
        if (xtest >= buf) & (xtest < xsz-buf) & (ytest >= buf) & (ytest < ysz-buf):
            phot = djs.djs_phot(ytest,xtest,ap,skyrad,image,skyrms=True,flerr=True,skyval=True,cbox=10)
            flux[i]   = phot["flux"]
            flerr[i]  = phot["fluxerr"]
            skyval[i] = phot["skyval"]
            skyrms[i] = phot["skyrms"]

    return jd,xval,yval,fwhm,aspect,snr,exptime,flux,flerr,skyval,skyrms



#----------------------------------------------------------------------#
# batch_phot:
#      calculate photometry on list of files with designated file
#      descriptors given ra and dec positions of stars
#----------------------------------------------------------------------#

def batch_phot(files,ras,decs,bias=None,dark=None,flat=None,outdir='./',
               datafile='photometry.pck',clobber=False,aperture=None,
               skyrad=np.array([20,30]),cr_reject=False):
    '''
    If aperture is left unset, optimal aperture will be used
    '''

    # Don't redo astrometry unless clobber keyword set
    if len(glob.glob(outdir+datafile)) == 1 and not clobber:
        print("Photometry file: "+datafile+" already exists!")
        print("Importing information...")
        info = pickle.load( open( outdir+datafile, "rb" ) )
        return info

    # File count
    fct = len(files)

    # Dimensions of photometry arrays
    vshape = (fct,len(ras))

    # Create vectors to feed final dictionary
    jd = np.empty(fct) ; x  = np.empty(fct) ; y  = np.empty(fct)
    fwhm = np.empty(fct) ; aspect = np.empty(fct)
    snr = np.empty(fct) ; exptime = np.empty(fct) ; apvec = np.empty(fct) 
    flux = np.empty(vshape) ; fluxerr = np.empty(vshape)
    sky = np.empty(vshape) ; skyrms = np.empty(vshape)
    peaks = np.empty(vshape)
    airmass = np.empty(fct)
    plate_scale = np.empty(fct)

    # Loop through files
    for i in np.arange(fct):
        print("Starting photometry of file: "+files[i].split('/')[-1])
        status = check_ast(files[i])
        #if status == 0:
        try:
            dict = do_phot(files[i],ras,decs,bias=bias,dark=dark,flat=flat,
                           aperture=aperture,skyrad=skyrad,peak=True)
            print "Measured flux = %.2f" % dict['flux'][0]
            jd[i]          = dict['jd']
            airmass[i]     = dict['secz']
            x[i]           = dict['xcen']
            y[i]           = dict['ycen']
            plate_scale[i] = dict['plate_scale']
            fwhm[i]        = dict['fwhm']
            aspect[i]      = dict['aspect']
            snr[i]         = dict['snrmax']
            exptime[i]     = dict['exptime']
            apvec[i]       = dict['aperture']
            flux[i,:]      = dict['flux']
            fluxerr[i,:]   = dict['flux_err']
            sky[i,:]       = dict['skyval']
            skyrms[i,:]    = dict['skyrms']
            peaks[i,:]     = dict['peak']
        
        except:
            print 'Photometry failure!!'

    # Make information dictionary
    info = {"jd":jd, "x":x, "y":y, "fwhm":fwhm, "aspect":aspect, "snr":snr,
            "exptime": exptime, "flux": flux, "fluxerr": fluxerr,
            "sky": sky, "skyrms": skyrms, "aperture": apvec,'airmass':airmass,
            "plate_scale": plate_scale, "peak":peaks}

    # Write out photometry data
    file = open(outdir+datafile, "w")
    pickle.dump(info, file)

    return info



#----------------------------------------------------------------------#
# light_curve:
#      produce light curve from list of files
#----------------------------------------------------------------------#

def lightcurve(info,title='Light Curve',
               bias=None,dark=None,order=3,write=True,
               outdir='./',datafile='photometry.pck',
               clobber=False,transit=False,userefs=False,
               flag=False):

    # Don't redo lightcurve unless clobber keyword set
    datafile = 'lightcurve.txt'
    if len(glob.glob(outdir+datafile)) == 1 and not clobber:
        print("Light curve file: "+datafile+" already exists!")
        print("Importing information...")
        jd, lc, err, exptime = np.loadtxt(outdir+datafile,unpack=True)
        return jd, lc, err, exptime

    flux    = info["flux"]
    exptime = info["exptime"]
    flerr   = info["flerr"]
    npts,nrefs = flux.shape
    ftarg  = flux[:,0]

# Use specified reference stars, else use reference stars that have data
# throughout dataset.
    n = 0
    fref = []
    if userefs != False:
        nrefs = len(userefs)
        for ref in userefs:
            ftest = flux[:,ref]
            inds, = np.where(ftest == 0)
            ftest[inds] = np.nan
#            if np.sum(ftest == 0) == 0:
            n += 1
            fref = np.append(fref,ftest,axis=0)
    else:
        for i in np.arange(nrefs-1)+1:
            ftest = flux[:,i]
            inds, = np.where(ftest == 0)
            ftest[inds] = np.nan
#            if np.sum(ftest == 0) == 0:
            n += 1
            fref = np.append(fref,ftest,axis=0)

    fref.shape = (n,npts)

# Create array of flux ratios
    fratios = fref
    for i in np.arange(n):
        fr = fref[i,:].flatten()
        fratios[i,:] = ftarg/fr

    test = np.nanmean(fratios,axis=0)

# Use beginning Julian day as the reference time
    jd = info["jd"]
    jdref = np.floor(np.min(jd))
    jdred = jd - jdref

# Plot up differential photometry
    nsub = str(np.int(ceil(np.sqrt(n))))
    fig = plt.figure(9)
    plt.clf()
    for i in np.arange(n):
        ax = fig.add_subplot(nsub,nsub,i+1)
        fplt = fratios[i,:]
        plt.plot(jdred,fplt/np.median(fplt),'o-',markersize=5,label='Ref '+str(i+1))
        plt.legend(loc=1)
    plt.draw()

    if not flag:
        flagval = raw_input("Would you like to flag the data?: ")
        if flagval == 'y': flag = True

    if flag:
        flagpts = []
        ax = fig.add_subplot(111)
        plt.clf()
        plt.plot(jdred,test,'o-')
        def flagclick(event):
            xflag = event.xdata
            yflag = event.ydata
            diff = np.sqrt((xflag-jdred)**2.+(yflag-test)**2)
            pt = np.argmin(diff)
            flagpts.append(pt)
            plt.scatter(jdred[pt],test[pt],marker='o',s=100,facecolor='red',edgecolor='red')
            plt.draw()

# Engage "baseclick"
        cid = fig.canvas.mpl_connect('button_press_event', flagclick)

# Stall here such that canvas can disconnect before further calculations
        print("Click on points to flag")
        print(raw_input("Press RETURN when finished"))

# Disengage "baseclick"

        fig.canvas.mpl_disconnect(cid)

        flagpts = np.unique(flagpts)

        fratios = np.delete(fratios,flagpts,axis=1)
        jd = np.delete(jd,flagpts)
        jdred = np.delete(jdred,flagpts)
        exptime = np.delete(exptime,flagpts)
        npts -= len(flagpts)

# If there is a full transit, give option to select "out of transit"
# region for baseline calculation.

    if transit:
        xvals=[]
        def baseclick(event):
            newx = event.xdata
            xvals.append(newx)
            print("------------------------------")
            print("xval = " , newx)
            plt.axvline(x=newx,linestyle='--',color='k')
            plt.draw()

# Engage "baseclick"
        cid = fig.canvas.mpl_connect('button_press_event', baseclick)

# Stall here such that canvas can disconnect before further calculations
        print("Click beginning and end of 2 baseline regions (4 selections)")
        print(raw_input("Press RETURN when finished"))

# Disengage "baseclick"
        fig.canvas.mpl_disconnect(cid)

# Check baseline region designations
        if len(xvals) > 4:
            print("Baseline regions overspecified!")
            print("Taking first 4 values")
            xvals = xvals[0:3]
        if len(xvals) < 4:
            print("Baseline regions under specified!")
            print("Returning...")
            return
        finds = np.where((jdred >= xvals[0]) & (jdred <= xvals[1]) |
                         (jdred >= xvals[2]) & (jdred <= xvals[3]))
    else:
        finds = np.arange(npts)


# Take mean time as zeropoint for polynomial fits to baseline regions
    ffinal = fratios*0
    weights = np.zeros(n)
    jdfit = jdred - np.mean(jdred)
    for i in np.arange(n):
        ffit = fratios[i,:]
        plt.clf()
        plt.subplot(2,1,1)
        plt.plot(jdfit,ffit,'.')
        plt.plot(jdfit[finds],ffit[finds],'.',color='r')
        plt.xlim([np.min(jdfit),np.max(jdfit)])
        fit = np.polyfit(jdfit[finds],ffit[finds],order)
        plt.title("Reference Star "+str(i+1))
        plt.ylabel("Flux Ratio")
        plt.xlabel("Time from Mean (d)")
        yvals = np.polyval(fit,jdfit[finds])
        rms = np.std(ffit[finds] - yvals)
        ymodel = np.polyval(fit,jdfit)
        plt.plot(jdfit,ymodel,'g-')
        plt.subplot(2,1,2)
        ffinal[i,:] = ffit/ymodel
        weights[i] = 1.0/(rms**2)
        plt.plot(jdred,ffit/ymodel,'.')
        plt.ylabel("Normalized Flux")
        plt.xlabel('JD - '+str(jdref)+' (days)')
        plt.xlim([np.min(jdred),np.max(jdred)])
        plt.annotate(r'$\sigma$ = %.4f' % rms, [0.85,0.4],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
        plt.draw()
        print(raw_input("Press RETURN to continue... \n"))


    lc = np.sum(ffinal*weights.reshape(n,1),axis=0)/np.sum(weights)

# Errors and exposure times
    err = np.zeros(len(lc)) + 1./(np.sqrt(np.sum(weights)))


    plt.figure(10)
    plt.clf()
    plt.plot(jdred,lc,'.')
    plt.xlim(np.min(jdred),np.max(jdred))
    plt.xlabel('JD - '+str(jdref)+' (days)')
    plt.ylabel('Normalized Flux')
    plt.title(title)
    plt.annotate(r'$\sigma$ = %.4f' % rb.std(lc[finds]), [0.8,0.8],horizontalalignment='right',
                 xycoords='figure fraction',fontsize='large')
    if write:
        plt.savefig(outdir+''.join(title.split())+".png",dpi=300)

# Write out final lightcurve
    outdata = np.zeros((len(lc),4))
    outdata[:,0] = jd
    outdata[:,1] = lc
    outdata[:,2] = err
    outdata[:,3] = exptime
    np.savetxt(outdir+"lightcurve.txt",outdata)

    return jd, lc, err, exptime


def get_info_path(verbose=False):
    '''
    Routine to check for the existence of the necessary directories of
    data to do throughput calculations
    '''
    a = sys.path
    b = None

    for dir in a:
        if dir[-11:] == '/photometry': b=dir

    if not b:
        b = '.'

    contents =  os.listdir(b)

    dpath = None
    stpath = None
    for c in contents:
        if c == 'PhotInfo':
            dpath = b+'/'+c
        if c == 'StellarTemplates':
            stpath = b+'/'+c

    if verbose:
        if not dpath:
            print 'PhotInfo directory not accessible'
        if not stpath:
            print 'StellarTemplates directory not available'

    return dpath, stpath


def seeing(info):

    pass



def Gaussian2D((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''
    Function to create a 2D Gaussian from input parameters.
    '''
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo)
                            + c*((y-yo)**2)))
    return g.ravel()



def moments(data):
    '''
    Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments
    '''

    xshape,yshape = data.shape
    med = np.nanmedian(data)
    datam = data - med
    total = np.nansum(datam)
    X, Y = np.indices(datam.shape)
    x_cen = np.nansum(X*datam)/total
    y_cen = np.nansum(Y*datam)/total
    x = np.float(xshape/2) if x_cen > xshape else x_cen
    y = np.float(xshape/2) if y_cen > yshape else y_cen
    col = datam[:, int(y)]
    width_x = np.sqrt(np.nansum(abs((np.arange(col.size)-y)**2*col))/np.nansum(col))
    row = datam[int(x), :]
    width_y = np.sqrt(np.nansum(abs((np.arange(row.size)-x)**2*row))/np.nansum(row))
    height = np.nanmax(datam)

    return height, x, y, width_x, width_y, med


def airmass(alt):
    #cos of zenith angle
    cos_zt = np.cos(np.radians(90.0 - alt))
    airmass = ((1.002432*cos_zt**2) + (0.148386*cos_zt) + 0.0096467 ) / \
               ( (cos_zt**3) + (0.149864*cos_zt**2) + (0.0102963*cos_zt) + 0.000303978)

    #airmass_approx = 1./cos_zt

    return airmass #,airmass_approx

    
def test_airmass_calcs(file):
    image,header = fits.getdata(file, 0, header=True)

    exptime = header['EXPTIME']

    jd = header["jd"] + (header["exptime"]/2.0)/(24.0*3600.0)
    obs.date = ephem.date(jd-2415020.0) # pyephem uses Dublin Julian Day (why!!??!?)
    star = ephem.FixedBody()

    # Should update airmass calculation from secz
    star._ra = np.radians(ra)
    star._dec = dec
    star.compute(obs)
    airmass1 = header['airmass']
    airmass2 = airmass(star.alt.degrees)

    return airmass1,airmass2


def par_angle(ha,dec,latitude=obs.lat):
    '''
    Returns the parallactic angle of a source given the hour angle, declination,
    and latitude of the observatory (default Thacher)

    PA is returned in degrees

    Stolen from T. Robishaw parangle
    http://www.cita.utoronto.ca/~tchang/gbt/procs/pangle/parangle.pro
    Thanks, Tim!
    '''
    har = ha*np.pi/180.0
    decr = dec*np.pi/180.0
    pa = -180/np.pi*\
         np.arctan2(-np.sin(har),
                 np.cos(decr)*np.tan(latitude)-np.sin(decr)*np.cos(har))
    return pa

def mech_pos(ha,dec,posang,latitiude=obs.lat,offset=181.0):
    '''
    Uses par_angle and position angle of an image to calculate 
    the mechanical angle of the rotator

    ha  = -1*(3+24.0/60.+41.31/3600.)*15.0
    dec = 48+1/60.0+43.2/3600
    pa  = 180-0.156
    offset = 181.0
    '''

    parangle = par_angle(ha,dec)
    mech = parangle + offset - posang
    
    if mech < 0:
        mech += 360.0

    return mech

def posang(header,verbose=False):
    CD11 = float(header['CD1_1'])
    CD12 = float(header['CD1_2'])
    CD21 = float(header['CD2_1'])
    CD22 = float(header['CD2_2'])

    ## This is my code to interpet the CD matrix in the WCS and determine the
    ## image orientation (position angle) and flip status.  I've used it and it
    ## seems to work, but there are some edge cases which are untested, so it
    ## might fail in those cases.
    ## Note: I'm using astropy units below, you can strip those out if you keep
    ## track of degrees and radians manually.
    if (abs(CD21) > abs(CD22)) and (CD21 >= 0): 
        North = "Right"
        positionAngle = 270.*u.deg + np.degrees(np.arctan(CD22/CD21))*u.deg
    elif (abs(CD21) > abs(CD22)) and (CD21 < 0):
        North = "Left"
        positionAngle = 90.*u.deg + np.degrees(np.arctan(CD22/CD21))*u.deg
    elif (abs(CD21) < abs(CD22)) and (CD22 >= 0):
        North = "Up"
        positionAngle = 0.*u.deg + np.degrees(np.arctan(CD21/CD22))*u.deg
    elif (abs(CD21) < abs(CD22)) and (CD22 < 0):
        North = "Down"
        positionAngle = 180.*u.deg + np.degrees(np.arctan(CD21/CD22))*u.deg
    if (abs(CD11) > abs(CD12)) and (CD11 > 0): East = "Right"
    if (abs(CD11) > abs(CD12)) and (CD11 < 0): East = "Left"
    if (abs(CD11) < abs(CD12)) and (CD12 > 0): East = "Up"
    if (abs(CD11) < abs(CD12)) and (CD12 < 0): East = "Down"
    if North == "Up" and East == "Left": imageFlipped = False
    if North == "Up" and East == "Right": imageFlipped = True
    if North == "Down" and East == "Left": imageFlipped = True
    if North == "Down" and East == "Right": imageFlipped = False
    if North == "Right" and East == "Up": imageFlipped = False
    if North == "Right" and East == "Down": imageFlipped = True
    if North == "Left" and East == "Up": imageFlipped = True
    if North == "Left" and East == "Down": imageFlipped = False

    pa = positionAngle.to(u.deg).value
    if North=='Up' and imageFlipped:
        pa = np.abs(pa-180)
    if verbose:
        print("Position angle of WCS is {0:.1f} degrees.".format(positionAngle.to(u.deg).value))
        print("Image orientation is North {0}, East {1}.".format(North, East))
        if imageFlipped:
            print("Image is mirrored.")

    return pa
