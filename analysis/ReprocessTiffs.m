% for re-writing tiffs from corrupted files!
%%i/o
clear
dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
files_in=cellstr(ls([dir_in, '*.tif']))

dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final2\';

dir_georef='F:\AboveDCSRasterManagement\CanadaAlbersTranslate\';
georef_files=cellstr(ls([dir_georef, '*.tif']))

%% loop
for i=149:length(files_in)
    tic
    fprintf('File %d:\t%s', i, files_in{i})
    path_in=[dir_in, files_in{i}];
    georef_path_in=[dir_georef, georef_files{i}];
%     [WC, R]=geotiffread(path_in);
    WC=imread(path_in);
    info=geotiffinfo(georef_path_in);
%     imshow(paths_out{i}, 'DisplayRange', [0 1]); %colormap parula
%     WC(WC==-99)=-1;
    figure(1); imagesc(WC, [0,1]); axis image
%     figure(2); imshow(paths{i})
    fprintf('\n')
    path_out=[dir_out, files_in{i}];
    try
        fprintf('Saving to directory: %s\n', dir_out)
        geotiffwrite(path_out, WC,...
            info.SpatialRef, 'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    catch
        disp('File not saved.')
    end
    toc
end