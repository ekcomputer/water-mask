## DESCRIPTION
This repository contains scripts for producing and analyzing open water classifications of
color-infrared imagery used in the following publication:

Kyzivat, E.D., et al, in prep

Data used for this paper can be found at the [Oak Ridge National Lab Distributed Active Archive Center (ORNL DAAC)](https://doi.org/10.3334/ORNLDAAC/1707)

The scripts are designed for AirSWOT color-infrared (CIR) images as produced by the July and August
2017 AirSWOT airborne sensor flights during the ongoing NASA Arctic-Boreal Vulnerability Experiment (ABoVE), 
an ongoing airborne and field-based campaign to study changes in Arctic-Boreal Alaska and Canada.

## Requirements

## FAQ

These scripts are designed to be run with custom file paths in the inputs and outputs- you will need to change what has been written in the following scripts: [BP_OBIA.m](#BP_OBIA.m), [BP_loadData.m](#BP_loadData.m), and [OBIA_BP_Fun.m](#OBIA_BP_Fun.m).  If running in test/debug mode, you will also need to edit: [BP_OBIA_Devel.m](#BP_OBIA_Devel.m).


## User parameters

The following parameters should be tuned by trial and error to acheive best results.  Default values (used for 2019 Remote Sensing paper) are provided as a good starting point.

Table 1: Explanation of parameters used for automated classification.  Values in brackets indicate parameters requiring manual selection for each image.  Note: DN is digital number.

| Parameter  | Description | Value |
|---      |---                |---       |
pArea | Pixel area in m2	| 1
minSize	| Min water region size (inclusive) in m2	| 40
satPercent |	Image enhancement after initial water index band math |	0.002
indexShrinkLim |	Max water index value (as a fraction of the global threshold) for erosion operation |	1.5
sz |	Target cluster size (pixels) |	100
minGrowSz |	Min number of clusters in region to allow region growing |	5
df |	?F, or expected flatness deviation as percent of max Euler number, used for binarization.  Has little effect when used for binarization, rather than multi-thresholding. |	20
cConn |	Step size for binarization (DN) |	[2, 3]
growMax |	Max number of region growing iterations (prevents endless growing) |	30
maxStd |	Max quantile of standard deviation for standard-deviation-based growing bounds. |	0.99
minAreaFact |	Number of times to multiply min cluster size (in meters) to determine number of extreme-valued pixels used for bounds for initial water determination (higher includes more extreme pixels). |	300
regionsLim | Max number of regions to allow growing.  Safety guard against endless growing. |	800
maxDilation | Max number of new clusters to be added to region during growing (prevents endless growing) |	500
bias |	Bias is added to the global threshold, so negative values make threshold lower. |	[-3, 0]
aConn |	Min threshold for binarization (units: water index as DN) |	[30, 230]
bConn |	Max threshold for binarization (units: water index as DN) |	[70, 250]
wp |	Sliding window size for binarization, expressed as percentage.  Be sure to make this small (i.e. 2) if the difference between bConn and aConn is small (i.e. 10) |	[2, 10]
windex |	Water index to use |	[NDWI, IR]
boundsLower |	Lower region growing bounds (expressed as multiple of the region’s standard deviation).  Higher values allow for more growing. |	[0.8, 2]
boundsUpper |	Upper region growing bounds (has very little effect) |	2.5
TLim |	Texture index cutoff. Lower values erode more heavily. | [5.3, 6.4]
NDWIWaterAmount |	Value of pixels above cutoff to show tile has water (units: water index).  Used to throw away tiles before classification if they don't have water. |	[0.04, 0.06] for NDWI, [0.50, 0.62] for IR
NDWILandAmount |	Value of pixels above cutoff to show tile has land (units: water index).  Used to throw away tiles before classification and mark as 100% water if they don't have land (rarely applies). |	[-1, -0.06] for NDWI, [-0.06, 0.48] for IR
TileSizeX |	Number of columns in processing tile |	[5,760, 16,000]
TileSizeY |	Number of rows in processing tile |	[1,600, 16,000]
RegionGrowing | Turn on for final classification.  Turn off for speed while testing. | [0 1]
Parallel | Turn on to use parallel processing.  Must install Matlab parallel processing toolbox first. | [0 1]

## Reference Material
We gratefully acknowledge the following papers and scripts for enabling the material here:
* [Image Graphs, version 1.0 (219 KB) by Steve Eddins](https://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs?focused=5570984&tab=example)
* [Power-law Distributions in Empirical Data](http://tuvalu.santafe.edu/~aaronc/powerlaws/)
* Campbell, J.B.; Wynne, R.H. Introduction to Remote Sensing; 2nd ed.; The Guilford Press: New York, 2011
* [Exact minimum bounding spheres/circles](https://www.mathworks.com/matlabcentral/fileexchange/48725-exact-minimum-bounding-spheres-circles)
* Rosin, P.L.; HervÃ¡s, J. Remote sensing image thresholding methods for determining landslide activity. Int. J. Remote Sens. 2005, 26, 1075â€“1092.
* [Tight subplot (matlab plotting tool)](https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
* [1. O’Gorman, L. Binarization and Multithresholding of Document Images Using Connectivity. CVGIP Graph. Model. Image Process. 1994, 56, 494–506.](dx.doi.org/10.1006/CGIP.1994.1044)

## Explanation of scripts
* [BP_OBIA.m](#BP_OBIA.m)
* [BP_OBIA_Devel.m](#BP_OBIA_Devel.m)
* [clip-to-polygon.sh](#clip-to-polygonsh)
* [BP_bigTiffWriterEK.m](#BP_bigTiffWriterEK.m)
* [BP_loadData.m](#BP_loadData.m)
* [GorminThreshold.m](#GorminThreshold.m)
* [OBIA_BP_Fun.m](#OBIA_BP_Fun.m)
* [RosinThreshold.m](#RosinThreshold.m)
* [SP_dil.m](#SP_dil.m)
* [SP_plot_raster.m](#SP_plot_raster.m)
* [fillLabelGaps.m](#fillLabelGaps.m)
* [growUntil.m](#growUntil.m)
* [imfillNaN.m](#imfillNaN.m)
* [mergeRegions_simple.m](#mergeRegions_simple.m)
* [optomizeConn.m](#optomizeConn.m)
* [regionFill.m](#regionFill.m)
* [sizeFilter.m](#sizeFilter.m)
* [waterindex.m](#waterindex.m)

# BP_OBIA.m

Script to apply block processing to water classification.  Had some bugs
with Matlab r2018a and AMD Radeon Pro Graphics card driver- causing a
memory error and the blue screen of death.  Bug is now possibly fixed by 
updating AMD driver.  Workaround at the time was to use opengl for
graphics operations rather than AMD driver.  BP stands for block
processing and OBIA stands for object-based image classification.

# BP_OBIA_Devel.m

This script calls OBIA_BP_Function in developer mode, meaning it
operatesz in single, smasll images without processing in blocks.

[BP_bigTiffWriterEK.m]

BIGTIFFWRITER - A basic image adapter to write Big TIFF files.
modified by EK for AirSWOT CIR images.

A simple ImageAdapter class implementing a Tiff writer ImageAdapter
object for use with BLOCKPROC. 
- Tile dimensions must be multiples of 16.
- Only uint8, RGB input image data is supported.

Based on "Working with Data in Unsupported Formats"
http://www.mathworks.com/help/toolbox/images/f7-12726.html#bse_q4y-1

# BP_loadData.m

This script plots input imiage and calls a script to compute a water
index such as the normalized-difference water index (NDWI).  It also
decides of the image contains all water or all land before enhancing and
rescaling the image to UINT8 format.  These decisions are contained in
the waterFlag first parameter (0=no water, 1= some water and land, 2= all
water).  The waterFlag second and third parameters are the median values
of the upper and lower n pixels of the image histogram, where n is a
user-supplied multiple (f.minAreaFact) of the smallest water body size.

# GorminThreshold.m

Function to find plateau pts of decaying exponential histogram, as
described in:
O’Gorman, L. Binarization and Multithresholding of 
document Images Using Connectivity. CVGIP Graph. Model. Image Process. 
56, 494–506 (1994).

Hist_counts is input histogram/PMF (can be offset horizontally)
plateaus is 1 x n vector of locations (usually just one) to binarize
image.

wp is sliding window size, expressed as percentage of
maximum image intensity.  Should be large to reduce noise, but not larger
than 'min' intensity diff bw levels.
df is deltaF, or expected flatness deviation as percent of max eul.
Doesn't change output unless I'm searching for peaks, rather than max.
Lower values make peaks more distinct, higher values combine peaks.
dyn_range is aprox dynamic range of image, (close to 256 for uint8) used
to compute width of sliding window.

# OBIA_BP_Fun.m

Main script to use OBIA and superpixels to classify open water extent
Output is final classified image.  Struct_in is block processing,
structure containing processing tile, log_dir gives location to write
processing stats. 
varargin should be the string 'local', which means to use region growing
to compute a local threshold for each water body.
includes entropy filter
Rewriting to include region growing/shrinking
Rewritten to detect SP on masked image

# RosinThreshold.m

best_idx = RosinThreshold(hist_img)

implementation of Rosin Thresholding, used for comparison with other image thresholding
techniques used for this project. 
Compute the Rosin threshold for an image
Takes histogram of an image filtered by any edge detector as as input
and return the index which corresponds to the threshold in histogram

REF: "Unimodal thresholding" by Paul L. Rosin (2001)

# SP_dil.m

dil_sps=SP_dil(g, SP_incl)
"Superpixel dilation": dilates a superpixel 'image' (graph) for the region containing labled
sps 'SP_incl' (vector), using graph g (of initial water SPs)

# SP_plot_raster.m

Lnew=SP_plot_raster(SP, L_all, {comparison, thresh}, 'complete')
Plotting utility to visualize connected components image stored in graph
form.  Output is a conversion back to raster (matrix) format.
optional arguments: [string] comparison is either 'lessthan' or 
'greaterthan'; last argument: 'noplot' doesn't plot
[double] thresh is max threshold for viewing non-index SP vectors
works for SP that includes zeros (mask vector) (filters them out)
imagesc for lists of superpixel SP (values are indices to L_all,
given SP label matrix L_all
returns Lnew, a label matrix of conglomerates of SP corr to water
only shows binary plots, colorized by SP index
warning: may change data type/class...

# fillLabelGaps.m

L=fillLabelGaps(L)
takes labeld matrix L (double, typical output of bwlabel) and if there are gaps in label indexes,
i.e. 1,2,3,6,7,..., it shifts everything else down to fill gaps.
warning: memory intensive!  Works as long as range(L(:)) aprox equal to
twice length(unique(L))

# growUntil.m
# imfillNaN.m
# mergeRegions_simple.m
# optomizeConn.m
# regionFill.m
# sizeFilter.m
# waterindex.m


```
  In development
```
