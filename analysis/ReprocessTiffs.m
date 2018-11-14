% for re-writing tiffs from corrupted mask (WM) files!
% uses georeff info from non-corrupt image files.
%%i/o
clear
% dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final3\';

files_in=cellstr(ls([dir_in, '*.tif']))

% dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final2\';
dir_out=dir_in;

dir_georef='F:\AboveDCSRasterManagement\CanadaAlbersTranslate\';
georef_files=cellstr(ls([dir_georef, '*.tif']))

    % make queue
Q=17:29;
%% loop
for i=Q
    tic
    fprintf('File %d:\t%s', i, files_in{i})
    path_in=[dir_in, files_in{i}];
    georef_name=['DCS', files_in{i}(3:end)]; % change name format to match image file
    georef_path_in=[dir_georef, georef_name]; % concat to directory
%     [WC, R]=geotiffread(path_in);
    WC=imread(path_in);
    info=geotiffinfo(georef_path_in); % read georeff info from corresponding georeferenced image file
%     imshow(paths_out{i}, 'DisplayRange', [0 1]); %colormap parula
%     WC(WC==-99)=-1;
    figure(1); imagesc(WC, [0,1]); axis image
%     figure(2); imshow(paths{i})
    fprintf('\n')
    path_out=[dir_out, files_in{i}];
        %write new file
    try
        fprintf('Saving to directory: %s\n', dir_out)
        geotiffwrite(path_out, WC,...
            info.SpatialRef, 'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    catch
        disp('File not saved.')
    end
    toc
end