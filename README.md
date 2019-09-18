## DESCRIPTION
This repository contains scripts for producing and analyzing open water classifications of
color-infrared imagery used in the following publication:

Kyzivat, E. D., L.C. Smith, S.W. Cooley, L.H. Pitcher, and C. Horvat et al. 2019. “A High-Resolution Airborne Color-Infrared Camera Water Mask for the NASA ABoVE Campaign.��? Remote Sensing 11. https://doi.org/https://doi.org/10.3390/rs11182163.

Data used for this paper can be found at the [Oak Ridge National Lab Distributed Active Archive Center (ORNL DAAC)](https://doi.org/10.3334/ORNLDAAC/1707)

The scripts are designed for AirSWOT color-infrared (CIR) images as produced by the July and August
2017 AirSWOT airborne sensor flights during the ongoing NASA Arctic-Boreal Vulnerability Experiment (ABoVE),
an airborne and field-based campaign to study changes in Arctic-Boreal Alaska and Canada.

These scripts turn 3-band color-infrared (CIR) images into binary masks denoting water, not water, and no data areas.  The analysis directory includes a variety of scripts for subsequent power-law distribution analysis, used in Kyzivat et al. 2019, as well as some additional procedures, such as area change analysis.

Thank you to [Chris Horvat](http://www.chrv.at/) for assistance with power-law analysis and plots.



## How to run

Set required environment vars (structure f), such as growing bounds, thresholding parameters, pathnames, etc., in OBIA_BP_Fun.m and env_vars.m.  Ensure your input data is in 3-band .tiff format or rewrite to load other file formats.  Run classifier using BP_OBIA.m (block processing on tiles from multiple images) or BP_OBIA_Devel.m (single images with no block processing), being sure to set user parameters at the beginning of these scripts, such as input and output paths.

## Requirements

- Matlab r2018a or similar version
- Toolboxes: Mapping Toolbox, Image Processing Toolbox, Statistics and Machine Learning Toolbox
- [Image Graphs, version 1.0 (219 KB) by Steve Eddins](https://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs?focused=5570984&tab=example) - be sure this toolbox is on your path
-If running analysis and plotting: [Tight Subplot](https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w)
- Optional: multi-core processor and Matlab Parallel Processing Toolbox.
- These scripts were run on a  Dell Precision Workstation with a 3.7 GHz Intel Xeon processor with 16 GB of RAM, resulting in a reasonably fast processing time of about four minutes for a 50-million-pixel processing tile.

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
boundsLower |	Lower region growing bounds (expressed as multiple of the region�s standard deviation).  Higher values allow for more growing. |	[0.8, 2]
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
* Rosin, P.L.; Hervas, J. Remote sensing image thresholding methods for determining landslide activity. Int. J. Remote Sens. 2005, 26, 1075-1092.
* [Tight subplot](https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w) (matlab plotting tool)
* Connectivity-based image thresholding:\
O'Gorman, L. Binarization and Multithresholding of Document Images Using Connectivity. CVGIP Graph. Model. Image Process. 1994, 56, 494-506. https://doi.org/10.1006/CGIP.1994.1044.
* Power-law analysis and plotting:\
Horvat, Christopher, Lettie Roach, Rachel Tilling, Cecilia Bitz, Baylor Fox-Kemper, Colin Guider, Kaitlin Hill, Andy Ridout, and Andrew Shepherd. 2019. “Estimating The Sea Ice Floe Size Distribution Using Satellite Altimetry: Theory, Climatology, and Model Comparison." The Cryosphere Discussions. https://doi.org/10.5194/tc-2019-134.

## Explanation of scripts
* [BP_OBIA.m](#bp_obiam)
* [BP_OBIA_Devel.m](#bp_obia_develm)
* [BP_bigTiffWriterEK.m](#bp-bigtiffwriterekm)
* [BP_loadData.m](#bp_loaddatam)
* [GorminThreshold.m](#gorminthresholdm)
* [OBIA_BP_Fun.m](#obia_bp_funm)
* [RosinThreshold.m](#rosinthresholdm)
* [SP_dil.m](#sp_dilm)
* [SP_plot_raster.m](#sp_plot_rasterm)
* [fillLabelGapsm](#filllabelgapsm)
* [growUntilm](#growuntilm)
* [imfillNaNm](#imfillnanm)
* [mergeRegions_simplem](#mergeregions_simplem)
* [optomizeConnm](#optomizeconnm)
* [regionFillm](#regionfillm)
* [sizeFilterm](#sizefilterm)
* [waterindexm](#waterindexm)
* [env_varsm](#env_varsm)

### BP_OBIA.m

Script to apply block processing to water classification.  Had some bugs
with Matlab r2018a and AMD Radeon Pro Graphics card driver- causing a
memory error and the blue screen of death.  Bug is now possibly fixed by
updating AMD driver.  Workaround at the time was to use opengl for
graphics operations rather than AMD driver.  BP stands for block
processing and OBIA stands for object-based image classification.

### BP_OBIA_Devel.m

This script calls OBIA_BP_Function in developer mode, meaning it operates in single, small images without processing in blocks.

### BP_bigTiffWriterEK.m

BIGTIFFWRITER - A basic image adapter to write Big TIFF files.
modified by EK for AirSWOT CIR images.

A simple ImageAdapter class implementing a Tiff writer ImageAdapter
object for use with BLOCKPROC.

- Tile dimensions must be multiples of 16.
- Only uint8, RGB input image data is supported.

Based on ["Working with Data in Unsupported Formats"](http://www.mathworks.com/help/toolbox/images/f7-12726.html#bse_q4y-1)

### BP_loadData.m

This script plots input image and calls a script to compute a water index such as the normalized-difference water index (NDWI).  It also decides of the image contains all water or all land before enhancing and rescaling the image to UINT8 format.  These decisions are contained in
the waterFlag first parameter (0=no water, 1= some water and land, 2= all
water).  The waterFlag second and third parameters are the median values
of the upper and lower n pixels of the image histogram, where n is a
user-supplied multiple (f.minAreaFact) of the smallest water body size.

### GorminThreshold.m

Function to find plateau pts of decaying exponential histogram, as
described in:
O�Gorman, L. Binarization and Multithresholding of
document Images Using Connectivity. CVGIP Graph. Model. Image Process.
56, 494�506 (1994).

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

### OBIA_BP_Fun.m

Main script to use OBIA and superpixels to classify open water extent
Output is final classified image.  Struct_in is block processing,
structure containing processing tile, log_dir gives location to write
processing stats.
varargin should be the string 'local', which means to use region growing
to compute a local threshold for each water body.
includes entropy filter
Rewriting to include region growing/shrinking
Rewritten to detect SP on masked image

### RosinThreshold.m

best_idx = RosinThreshold(hist_img)

implementation of Rosin Thresholding, used for comparison with other image thresholding
techniques used for this project.
Compute the Rosin threshold for an image
Takes histogram of an image filtered by any edge detector as as input
and return the index which corresponds to the threshold in histogram

REF: "Unimodal thresholding" by Paul L. Rosin (2001)
<a name="sp-dil">
### SP_dil.m
</a>

dil_sps=SP_dil(g, SP_incl)
"Superpixel dilation": dilates a superpixel 'image' (graph) for the region containing labled
sps 'SP_incl' (vector), using graph g (of initial water SPs)

### SP_plot_raster.m

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

### fillLabelGaps.m

L=fillLabelGaps(L)
takes labeld matrix L (double, typical output of bwlabel) and if there are gaps in label indexes,
i.e. 1,2,3,6,7,..., it shifts everything else down to fill gaps.
warning: memory intensive!  Works as long as range(L(:)) aprox equal to
twice length(unique(L))

### growUntil.m

complete_region=growUntil(g, spIncl, ~, sp_mean, ~, sp_std, sp_rcount, ~)

V1 uses dynamic programming
dilates a superpixel 'image' (graph), g, for the region containing
labled sps 'SP_incl' (vector), using graph g (of initial water SPs) until
condition is met
lim are limits of growing as multiple of the standard deviation of the
region (sp_std), as given by the global var f.bounds (ex: .9, 1.1).
note spIncl gives indexes to sp_mean;
add variance, entropy, or texture image
calls SP_dil() and fastSetdiff()
sp_rcount is vector giving sizes of every SP, in pixels.  Function
returns a new graph containg SP indexes for the new, dilated regoion.

### imfillNaN.m

Function to fill in NaN values surrounded by foreground
in classified, binary image
mask is binary mask with 1= no value (NaN value)

### mergeRegions_simple.m

Function to aggregate SPs (superpixels) bnased on an a priori mask, in
order to reduce the total number of SPs, and thus the size of the
dataset.
merges regions in L_all (label matrix) based on regions in bw (binary
image).  Returns simplified SP image (outputImage), vector SP_out,
SP_rcount (which is vector of number of combined SPs), and
updated label matrix (L_all_out)
varargin can be 'mean' (default), 'std' or 'count'
calls fillLabelGaps

### optomizeConn.m

[bw, loc]= optomizeConn(gray, ~, NoValues, bias)
Binaroization function baased on O'Gorman's (1994)
connectivity-preserving algoprithm for document image scanning.
bias moves thresh in the direction of bias (negative values move thresh
down)
Revised to measure  connectivity with greycomatrix no. of water regions
connected
Gray is grayscale image and should be superpixilated and uint8.  NoValues
is image mask, where 1 corresponds to no data.
Function decreases binary threshold until connectivity is maximized
output binary classified image

### regionFill.m

Region growing and shrinking algorithm.
shrinking is done using global thresh
region shrinking/growing on superpixels
sp_text is entropy image of p[rocessing tile, uased for image erosion.
cir_index is input image/processing tile, giving gray levels on which to base
growing/shrinking.
L_all is label matrix of SPs, bw is optimal connectivity mask
outputImage=ndwi or other index (uint8)
returns regiond, a list of SP corr to classified water, and
Lnew, a label matrix of conglomerates of SP corr to water

function calls: edgs2adjList, SP_dil, MergeRegionsSimple,
AdjacentRegionsGraph, SP_plot_raster, and growUntil.

### sizeFilter.m

bw_out=sizeFilter(bw, minSize)
removes regions in binary image below minSize (inclusive)
uses 8-connectivity
output is also binary

### waterindex.m

[cir_index]=waterindex(cir, waterIndex, ~)
Function to return various banmd ratios and indexes from an input image.
For this paper, the Normalized Diffewrence Water Ind3ex (NDWI) and
infrared (IR, band 1) was used.
cir is 3-band int image
cir_index is single precision
bw is binary threshold with certain sensitivity

### env_vars.m

environment variables
specify most recent versions of files
This file should run at start of any analysis session


```
  In development
```
