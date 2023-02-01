#From: https://docs.astropy.org/en/stable/wcs/wcstools.html#matplotlib-plots-with-correct-wcs-projection
from matplotlib import pyplot as plt
import numpy as np
from sys import exit
import imageio

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


frames = []
annFrames = []
##Print config understanding to user
for num in range(numFiles):
    fitsFileBase = fitsFileNames[num].split('.')[0]
    fileName = "./img/{}_{}.png".format('_'.join(title.split()), fitsFileBase)
    frames.append(fileName) 
    print("   {:d} {}".format(num, fileName))

print("These will be saved as ./img/{}.gif".format('_'.join(title.split())))

for num in range(numFiles):
    fitsFileBase = fitsFileNames[num].split('.')[0]
    annFileName = "./img/{}_annotated_{}.png".format('_'.join(title.split()), fitsFileBase)
    annFrames.append(annFileName) 
    print("   {:d} {}".format(num, annFileName))

print("These will be saved as ./img/{}_annotated.gif".format('_'.join(title.split())))

imgs = []
annImgs = []

for imgName in frames:    
    imgs.append(imageio.v2.imread(imgName))

for imgName in annFrames:
    annImgs.append(imageio.v2.imread(imgName))

imageio.mimsave("./img/{}.gif".format('_'.join(title.split())), imgs, duration=0.25, loop=0)
imageio.mimsave("./img/{}_annotated.gif".format('_'.join(title.split())), annImgs, duration=0.25, loop=0)
