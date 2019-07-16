% uses averaging windows to output stats and uses vector mask

clear

msk_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\shp\PixelMask_join_diss_simp.shp';
distrib_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\shp\distrib.shp';
[msk,msk_att]=shaperead(msk_in);

[distrib, distrib_att]=shaperead(distrib_in);


