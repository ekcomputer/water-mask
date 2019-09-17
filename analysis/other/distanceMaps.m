% script to calculate distance maps from water masks

dir_in_dist='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
% dir_in_dist='I:\Final\';

% dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\labeled\';
dir_out_dist='F:\AboveDCSRasterManagement\dist\';
dir_out_dist='D:\GoogleDrive\ABoVE top level folder\AirSWOT_CIR\dist\';

% dir_out_dist='I:\dist\';

files=cellstr(ls([dir_in_dist, '*.tif']));

%% loop
for i=85:length(files)
    pth_in=[dir_in_dist, files{i}];
    pth_out_dist=[dir_out_dist, 'D',files{i}(3:end)];
    gtinfo=geotiffinfo(pth_in);
    [WC, R]=geotiffread(pth_in);
    D=uint16(bwdist(WC));
    D(WC==-99)=65535; % re-mask nodata region
    geotiffwrite(pth_out_dist, D,...
        R, 'GeoKeyDirectoryTag',gtinfo.GeoTIFFTags.GeoKeyDirectoryTag);
    fprintf('Saving: %d: %s\n', i, files{i})
end