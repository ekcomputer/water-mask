## DESCRIPTION
This repository contains scripts for producing and analyzing open water classifications of
color-infrared imagery used in the following publication:

Kyzivat, E.D., et al, in prep

Data used for this paper can be found at the [Oak Ridge National Lab Distributed Active Archive Center (ORNL DAAC)](https://doi.org/10.3334/ORNLDAAC/1707)

## Requirements
You must install [**GDAL**](http://gdal.org/) and be able to access its utilities from the command line.

##FAQ

These scripts are designed to be run with custom file paths in the inputs and outputs- you will need to change what has been written in the following scripts: [BP_OBIA.m](#BP_OBIA.m), [BP_loadData.m](#BP_loadData.m), and [OBIA_BP_Fun.m](#OBIA_BP_Fun.m).  If running in test/debug mode, you will also need to edit: [BP_OBIA_Devel.m](#BP_OBIA_Devel.m).


## Reference Material
We gratefully acknowledge the following papers and scripts for enabling the material here:
* [Image Graphs, version 1.0 (219 KB) by Steve Eddins](https://www.mathworks.com/matlabcentral/fileexchange/53614-image-graphs?focused=5570984&tab=example)
* [Power-law Distributions in Empirical Data](http://tuvalu.santafe.edu/~aaronc/powerlaws/)
* Campbell, J.B.; Wynne, R.H. Introduction to Remote Sensing; 2nd ed.; The Guilford Press: New York, 2011
* [Exact minimum bounding spheres/circles](https://www.mathworks.com/matlabcentral/fileexchange/48725-exact-minimum-bounding-spheres-circles)
* Rosin, P.L.; Hervás, J. Remote sensing image thresholding methods for determining landslide activity. Int. J. Remote Sens. 2005, 26, 1075–1092.


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


```
  In development

```
TODO: Rosin RosinThreshold
