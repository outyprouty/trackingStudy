0. Needed packages: astropy, photutils, matplotlib, numpy, tabulate, imageio


1. ./config file
The config file is a flat ASCII file with following format.

```
<N> #integer number of fits-wcs pairs to follow
<FITS> <WCS> #space-delimited fits-wcs pair, relative to their enclosing directories with extensions
<M> #integer number representing the WCS file to use for all plotting
<STR> #title of images
<DAO|IRAF> #star finder algorithm (IRAF is slower, but can be used to report FWHMs reliably)
```

2. runTrackingStudy.py
RUN: `python runTrackingStudy.py`
This file reads the config file and reiterates the config to the user.
Then it uses astropy.io.fits.open to open the specified FITS HDUs along with their associated WCS files (more on how I got the WCS below).
Then it grabs the WCS to compare all the other HDU data to.
It makes plots of each HDU data segment w.r.t. the selected WCS; this is the first set of images.
It uses either DAO StarFinder (https://photutils.readthedocs.io/en/stable/api/photutils.detection.DAOStarFinder.html) OR
   IRAF Starfinder (https://photutils.readthedocs.io/en/stable/api/photutils.detection.IRAFStarFinder.html#photutils.detection.IRAFStarFinder)
   both are provided by photutils.detection package.
The found stars are queried for their image-system pixel location. These locations are annotated in the second set of images.
The pixel locations are translated via the selected WCS to that WCS system and associated with their image. I chose to use FK5. More on that below.
A list of non-symmetric image pairs are generated.
Then I report the differences in the WCS frame and image-system frame along with the averages for each pair considered.

Note: If the IRAF star finder is used, I also report the FWHM of the found sources.

3. img2gif.py
RUN: `python img2gif.py`
This code reads the same config file and uses the information to find the generated images.
It then generates two animations (GIF) from the images generated.

4. Below
WCS generation
I used an online tool "astrometry.net" to recover the best-guessed WCS for each frame.
Specifically, I used the API the developers make available for cloning.
To try for yourself: https://nova.astrometry.net/api_help AND/OR https://github.com/dstndstn/astrometry.net/tree/main/net/client



Fifth Fundamental Catalogue (FK5)
FK5 is part of the "Catalogue of Fundamental Stars" which provides a series of six astrometric catalogues of high precision positional data for a small selection of stars to define a celestial reference frame. J2000 refers to the instant of 12 PM (midday) on 1st January, 2000. FK5 was published in 1991 and added 3,117 new stars.
