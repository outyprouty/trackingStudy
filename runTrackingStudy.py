#From: https://docs.astropy.org/en/stable/wcs/wcstools.html#matplotlib-plots-with-correct-wcs-projection
import warnings
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.utils.data import get_pkg_data_filename
from datetime import datetime

from astropy.stats import sigma_clipped_stats

from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture
import numpy as np

from sys import argv, exit

firstNum = int(argv[1])
lastNum = int(argv[2])

fitsFileNames = []
wcsFileNames = []

# ./config is flat ASCII file following format in README
opts = [f.strip().split() for f in open("config", 'r').readlines()]

##begin issection of config file
numFiles = int(opts[0][0])
for pair_num in range(numFiles):
	fitsFileNames.append(opts[pair_num+1][0])
	wcsFileNames.append(opts[pair_num+1][1])

WCS_IDX = int(opts[numFiles+1][0])
WCS_FILENAME = wcsFileNames[WCS_IDX]
title = " ".join(opts[numFiles+2])

##end dissection of config file

##Print config understanding to user
print("Running Tracking Study")
print(" Found {:d} FITS files to consider".format(numFiles))
for num in range(numFiles):
	print("   {:d} {}".format(num, fitsFileNames[num]))

print(" Will compare all above relative to WCS specified by: {}".format(WCS_FILENAME))


# Open WCS
hdu_wcs = fits.open("./wcs/"+WCS_FILENAME)[0]
with warnings.catch_warnings():
    # Ignore a warning on using DATE-OBS in place of MJD-OBS
	warnings.filterwarnings('ignore')
	DATE_WCS = fits.open("./data/"+fitsFileNames[WCS_IDX])[0].header["DATE-OBS"]


RAs,DECs = [], []
X,Y = [],[]

#Cycle through each fits-wcs pair
for n in range(numFiles):
    ## Open FITS
    fitsFileBase = fitsFileNames[n].split('.')[0]
    hdu = fits.open("./data/"+fitsFileBase+".fits")[0]

    #Set WCS
    with warnings.catch_warnings():
        # Ignore a warning on using DATE-OBS in place of MJD-OBS
        warnings.filterwarnings('ignore')
        wcs = WCS(hdu_wcs.header)

    ## Constrain data to allow for edge buffer and ignore datestamp from SharpCap
    data = hdu.data[50:-50, 50:-50]

    #Grab some scene stats
    mean, median, std = sigma_clipped_stats(data, sigma=3.0) 

    ## Setup DAO Source Finder Alg
    daofind = DAOStarFinder(fwhm=17.0, threshold=median, brightest=10)

    ## Find sources 
    sources = daofind(data)


    ## Grab positions
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    ## Plot image in WCS Frame
    fig = plt.figure(figsize=(16,10))
    ax = plt.subplot(projection=wcs)
    ax.imshow(data, origin='lower', cmap="Greys")

    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')
    ax.coords.grid(True, color='white', ls='solid')

    ax.set_title("%s\n"%title+hdu.header["DATE-OBS"] +" wrt "+DATE_WCS)
    print("Image with WCS overlay saved to {}".format("./img/image_%s.png"%(fitsFileBase)))
    plt.savefig("./img/image_%s.png"%(fitsFileBase))
    plt.close()
        
    ## Plot image in WCS Frame with labelled sources 
    fig = plt.figure(figsize=(16,10))
    ax = plt.subplot(projection=wcs)
    ax.imshow(data, origin='lower', cmap="Greys")

    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')
    ax.coords.grid(True, color='white', ls='solid')

    ax.set_title("%s\n"%title+hdu.header["DATE-OBS"] +" wrt "+DATE_WCS)

    tmpRA, tmpDEC = [], []
    tmpX, tmpY = [], []
    for i in range(len(positions)):
        x, y = positions[i]
        ax.scatter(x,y, s=50, c='red', marker='x')
        ax.text(x+10,y+10, "%d"%i)
        coords = wcs.pixel_to_world(x,y)
        tmpX.append(x)
        tmpY.append(y)
        tmpRA.append(coords.ra.deg)
        tmpDEC.append(coords.dec.deg)
    RAs.append(tmpRA)
    DECs.append(tmpDEC)
    X.append(tmpX)
    Y.append(tmpY)

    print("Image with WCS overlay and annotations saved to {}".format("./img/image_annotated_%s.png"%(fitsFileBase)))
    plt.savefig("./img/image_annotated_%s.png"%(fitsFileBase))
    plt.close()



#'2023-01-25T01:10:21.2185211'
#Removed microseconds due to strptime not having microsecond support for more than 7 digits
lastDateStr = fits.open("./data/"+fitsFileNames[lastNum])[0].header["DATE-OBS"].split('.')[0]
firstDateStr = fits.open("./data/"+fitsFileNames[firstNum])[0].header["DATE-OBS"].split('.')[0]
lastDate = datetime.strptime(lastDateStr, "%Y-%m-%dT%H:%M:%S")
firstDate = datetime.strptime(firstDateStr, "%Y-%m-%dT%H:%M:%S")

totSec = (lastDate-firstDate).total_seconds()


print()
print("Begin Tracking Analysis")
print("  Below are the differences between Capture {:d} and Capture {:d}".format(lastNum, firstNum))
print("  Capture {:d}   {}".format(lastNum, lastDate))
print("  Capture {:d}   {}".format(firstNum, firstDate))
print("  The total time between the two dates is: {:f}s".format(totSec))
print()

RAs = np.array(RAs)
DECs = np.array(DECs)
X = np.array(X)
Y = np.array(Y)

diffR = ((RAs[2]-RAs[0])/totSec)*3600.

diffD = ((DECs[2]-DECs[0])/totSec)*3600

diffX = (X[2]-X[0])/totSec

diffY = (Y[2]-Y[0])/totSec


from tabulate import tabulate
table = np.array([range(10), diffR, diffD, diffX, diffY]).T
print(tabulate(table, headers=["Source ID", "dRA[''/s]", "dDEC[''/s]", "dX[p/s]", "dY[p/s]"],floatfmt=('.0f', '.6f', '.6f', '.6f', '.6f')))
print()
print(tabulate([["Average",np.average(diffR), np.average(diffD), np.average(diffX), np.average(diffY)]], headers=["         ", "dRA[''/s]", "dDEC[''/s]", "dX[p/s]", "dY[p/s]"], floatfmt=(None, '.6f', '.6f', '.6f', '.6f')))

print()
pixScale = float(wcs.proj_plane_pixel_scales()[0].to_value(unit='deg')*3600)
print("Calculated Pixel Size: {:0.4f} [''/pix]".format(pixScale))

#How to incorperate FWHM error aka seeing error?
#How to reconcile different measures of pixel dimensions?
