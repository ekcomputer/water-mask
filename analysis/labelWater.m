% script to produce numbered region or water maps with each body numbered

dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
% dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\labeled\';
dir_out='F:\AboveDCSRasterManagement\labeled\';
files=cellstr(ls([dir_in, '*.tif']));
labelDist=0;
if labelDist
    dir_out_dist='F:\AboveDCSRasterManagement\dist\';
end
%% loop
for i=257:length(files)
    pth_in=[dir_in, files{i}];
    pth_out=[dir_out, 'L',files{i}(3:end)];
    gtinfo=geotiffinfo(pth_in);
    [WC, R]=geotiffread(pth_in);
    L=uint16(bwlabel(WC));
    L(WC==-99)=0; % re-mask nodata region
    geotiffwrite(pth_out, L,...
        R, 'GeoKeyDirectoryTag',gtinfo.GeoTIFFTags.GeoKeyDirectoryTag);
    fprintf('Saving: %d: %s\n', i, files{i})
    
    if labelDist
        pth_out_dist=[dir_out_dist, 'D',files{i}(3:end)];
        D=uint16(bwdist(WC));
        D(WC==-99)=65535; % re-mask nodata region
        geotiffwrite(pth_out_dist, D,...
            R, 'GeoKeyDirectoryTag',gtinfo.GeoTIFFTags.GeoKeyDirectoryTag);
        fprintf('\tSaving distance map\n')
    end
end