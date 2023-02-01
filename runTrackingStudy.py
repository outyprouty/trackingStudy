#From: https://docs.astropy.org/en/stable/wcs/wcstools.html#matplotlib-plots-with-correct-wcs-projection
import warnings
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
from tabulate import tabulate
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder,IRAFStarFinder
import numpy as np
from sys import exit


# ./config is flat ASCII file following format in README
opts = [f.strip().split() for f in open("config", 'r').readlines()]

##begin issection of config file
fitsFileNames = []
wcsFileNames = []

numFiles = int(opts[0][0])
for pair_num in range(numFiles):
	fitsFileNames.append(opts[pair_num+1][0])
	wcsFileNames.append(opts[pair_num+1][1])

WCS_IDX = int(opts[numFiles+1][0])
WCS_FILENAME = wcsFileNames[WCS_IDX]
title = " ".join(opts[numFiles+2])
finder = " ".join(opts[numFiles+3])

##end dissection of config file




##Print config understanding to user
print("Running Tracking Study: {}".format(title))
print(" Found {:d} FITS files to consider".format(numFiles))
for num in range(numFiles):
	print("   {:d} {}".format(num, fitsFileNames[num]))
print( "Will use {} star finder algorithm to locate sources.".format(finder))
print(" Will compare all above relative to WCS specified by: {}".format(WCS_FILENAME))




# Open WCS
hdu_wcs = fits.open("./wcs/"+WCS_FILENAME)[0]

#Set WCS
#grab WCS observation date -- WCS file doesn't follow standarf FITS format, so we ignore an error here
with warnings.catch_warnings():
    # Ignore a warning on using DATE-OBS in place of MJD-OBS
    warnings.filterwarnings('ignore')
    wcs = WCS(hdu_wcs.header)
    DATE_WCS = fits.open("./data/"+fitsFileNames[WCS_IDX])[0].header["DATE-OBS"]


#Setup empty arrs to capture all angular and pixel locations found 
RAs,DECs = [], []
X,Y = [],[]

#Cycle through each FITS file, grab the image data, find the stars, log their locations, and make some plots
for n in range(numFiles):
    ## Open FITS
    fitsFileBase = fitsFileNames[n].split('.')[0]
    hdu = fits.open("./data/"+fitsFileBase+".fits")[0]

    ## Constrain data to allow for edge buffer and ignore datestamp from SharpCap
    data = hdu.data[50:-50, 50:-50]

    #Grab some scene stats
    mean, median, std = sigma_clipped_stats(data, sigma=3.0) 

    ## Setup DAO Source Finder Alg
    #Use either DAO or IRAF Star Finder from photutils
    # https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html
    # https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html
    # I naively assume positional arguments are treater the same way by both algs
    
    if finder == "DAO": 
        starFind = DAOStarFinder(fwhm=17.0, threshold=median, brightest=10)
    elif finder == "IRAF":
        starFind = IRAFStarFinder(fwhm=17.0, threshold=median, brightest=10)    
    else:
        print("This is embarassing ... No finder found")
        exit(0)

    ## Find sources, store results in sources dictionary
    sources = starFind(data)

    ## Extract just the positions
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    ## Plot image in WCS Frame
    fig = plt.figure(figsize=(16,10))
    ax = plt.subplot(projection=wcs)
    ax.imshow(data, origin='lower', cmap="Greys")

    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')
    ax.coords.grid(True, color='white', ls='solid')

    ax.set_title("%s\n"%title+hdu.header["DATE-OBS"] +" wrt "+DATE_WCS)
    fileName = "./img/{}_{}.png".format('_'.join(title.split()), fitsFileBase)
    print("Image with WCS overlay saved to {}".format(fileName))
    plt.savefig(fileName)
    plt.close()
        
    ## Plot image in WCS Frame with labelled sources 
    fig = plt.figure(figsize=(16,10))
    ax = plt.subplot(projection=wcs)
    ax.imshow(data, origin='lower', cmap="Greys")

    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')
    ax.coords.grid(True, color='white', ls='solid')

    ax.set_title("%s\n"%title+hdu.header["DATE-OBS"] +" wrt "+DATE_WCS)

    #add annotations to this second figure and collect positional information
    tmpRA, tmpDEC = [], []
    tmpX, tmpY = [], []
    for i in range(len(positions)):
        x, y = positions[i]
        #adding annotations
        ax.scatter(x,y, s=50, c='red', marker='x')
        ax.text(x+10,y+10, "%d"%i)
        coords = wcs.pixel_to_world(x,y)
        
        #Store pixel locations
        tmpX.append(x)
        tmpY.append(y)
        
        #USe WCS to convert to RA and DEC in degress, store these
        tmpRA.append(coords.ra.deg)
        tmpDEC.append(coords.dec.deg)

    #Add these temp arrs to globally defined arrs for later analysis
    RAs.append(tmpRA)
    DECs.append(tmpDEC)
    X.append(tmpX)
    Y.append(tmpY)

    fileName = "./img/{}_annotated_{}.png".format('_'.join(title.split()), fitsFileBase)
    print("Image with WCS overlay and annotation saved to {}".format(fileName))
    plt.savefig(fileName)
    plt.close()


#print pixel scale for reference
pixScale = float(wcs.proj_plane_pixel_scales()[0].to_value(unit='deg')*3600)
print()
print("Calculated Pixel Size: {:0.4f} [''/pix]".format(pixScale))
print()

print("Begin Tracking Analysis")
print()


#Generate list of pairs
idxs = range(numFiles)
pairs = [(i, j) for idx, i in enumerate(idxs) for j in idxs[idx + 1:]]

for p in pairs:
    firstNum, lastNum = p

    #Nominal date example: '2023-01-25T01:10:21.2185211'
    #Removed microseconds due to strptime not having microsecond support for more than 7 digits
    lastDateStr = fits.open("./data/"+fitsFileNames[lastNum])[0].header["DATE-OBS"].split('.')[0]
    firstDateStr = fits.open("./data/"+fitsFileNames[firstNum])[0].header["DATE-OBS"].split('.')[0]
    lastDate = datetime.strptime(lastDateStr, "%Y-%m-%dT%H:%M:%S")
    firstDate = datetime.strptime(firstDateStr, "%Y-%m-%dT%H:%M:%S")

    #Calculate total exposure time in seconds
    totSec = (lastDate-firstDate).total_seconds()


    print()
    print("  Below are the differences between Capture {:d} and Capture {:d}".format(lastNum, firstNum))
    print("  Capture {:d}   {}".format(lastNum, lastDate))
    print("  Capture {:d}   {}".format(firstNum, firstDate))
    print("  The total time between the two dates is: {:f}s".format(totSec))
    print()

    RAs = np.array(RAs)
    DECs = np.array(DECs)
    X = np.array(X)
    Y = np.array(Y)

    diffR = ((RAs[lastNum]-RAs[firstNum])/totSec)*3600.

    diffD = ((DECs[lastNum]-DECs[firstNum])/totSec)*3600

    diffX = (X[lastNum]-X[firstNum])/totSec

    diffY = (Y[lastNum]-Y[firstNum])/totSec

    table = np.array([range(len(RAs[0])), diffR, diffD, diffX, diffY]).T
    print(tabulate(table, headers=["Source ID", "dRA[''/s]", "dDEC[''/s]", "dX[p/s]", "dY[p/s]"],floatfmt=('.0f', '.6f', '.6f', '.6f', '.6f')))
    print()
    print(tabulate([["Average",np.average(diffR), np.average(diffD), np.average(diffX), np.average(diffY)]], headers=["         ", "dRA[''/s]", "dDEC[''/s]", "dX[p/s]", "dY[p/s]"], floatfmt=(None, '.6f', '.6f', '.6f', '.6f')))

    print()

    #How to incorperate FWHM error aka seeing error?
    #How to reconcile different measures of pixel dimensions?
