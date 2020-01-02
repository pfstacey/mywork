"""
@author: Yousef
"""
import thacherphot as tp
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u

'''
How to apply this script to any night
1. change flat,bias,dark paths
2. choose reference stars (tp.choose_refs() or manually through fits file viewer)
3. change graph labels for target name, photometric band, date, etc.
(4). set clobber to False when photometry already done (.pck file will have been made)
'''
#------------------------------------
#get files and photometric data
#------------------------------------

# Calibrations
files, fct = tp.get_files(dir='/Users/Yousef/Desktop/Astronomy/TransitData/2017Mar31/',prefix='WASP',tag='bias',suffix='.fit')
bias = tp.master_bias(files,verbose=False)
files,fct = tp.get_files(dir='/Users/Yousef/Desktop/Astronomy/TransitData/2017Apr10/',prefix='HAT',tag='flat',suffix='.fit')
flat = tp.master_flat(files,bias=bias)

# Check calibrated source file
files,fct = tp.get_files(dir='/Users/Yousef/Desktop/Astronomy/TransitData/2017Apr10/',prefix='HAT',tag='solved',suffix='.fits')

# Choose reference stars
plt.ion()
ras = [2.093894999999999698e+02, 2.093196788441016167e+02,2.092947705668225922e+02,2.093539466884411695e+02,2.093448754413380755e+02, 2.092702946383417100e+02, 2.095340387352351286e+02, 2.092523347142034993e+02] 
decs = [4.349352777777777845e+01,4.344867495658954226e+01,4.341084673274470163e+01,4.352678857733513240e+01,4.359690982054149089e+01,4.337773306989716815e+01,4.338057231726356378e+01,4.365228095082844817e+01]

# Do photometry on all files 263
phot = tp.batch_phot(files,ras,decs,bias=bias,flat=flat,skyrad=[20,25],clobber=False)

#-----------------------------------------
#reduce photometric data
#-----------------------------------------

numvals = len(phot['flux'][0])
numrefs = len(phot['flux'][0]) -1
numimages = len(phot['flux'])

#extract fluxes from dictionary
def get_fluxes():
    target = []
    for k in range(numimages):
        target.append(phot['flux'][k][0])
    
    refdata = []
    for refstar in range(1, numvals):
        temp = []
        for k in range(numimages):
            temp.append(phot['flux'][k][refstar])
        refdata.append(temp)
        
    return target,refdata

target,refdata = get_fluxes()

#get flux ratios
def get_fluxratios(target,refdata):
    fluxratios = []
    for refstar in range(numrefs):
        temp = [float(c)/t for c,t in zip(target, refdata[refstar])]
        fluxratios.append(temp)
        
    return fluxratios
    
fluxratios = get_fluxratios(target,refdata)

#normalize fluxes
def normalize_fluxes(fluxratios):
    normfluxes = []
    for refstar in range(numrefs):
        median = np.median(fluxratios[refstar][0:60])
        normfluxes.append(fluxratios[refstar]/median)
    
    return normfluxes

normfluxes = normalize_fluxes(fluxratios)

#average normalized fluxes    
total = 0
for k in range(numrefs):
    total += normfluxes[k]
    
finalflux = total/numrefs

plt.plot(phot['jd'],finalflux, 'o',markersize=3.5)
plt.title('HAT-p 12 b Transit (2017Apr10,band: V)')
plt.ylabel('Flux')
plt.xlabel('Julian Date')
