% todo: add geotiffwrite at end insted of .tfw; change bigtiff to tiff
% writer; edit optomizeConn_16 for case of no peak in MasterMetric and
% other scaling issues (incl no bimodal hist...)
% Script to apply block processing to DCS images- for water classification
% modified from DCS_images_4.m

% File queue
clear all; clc
tic
% set(0,'DefaultFigureVisible','off')
% dir_in='D:\ArcGIS\FromMatlab\ClipSquares\';
dir_in='D:\GoogleDrive\ABoVE top level folder\AirSWOT_CIR\DAAC_Preview\DCS_all\';
dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\';
logfile='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\logs\log.txt';
fid=fopen(logfile, 'a')
fprintf(fid, '----------------------\n')
fclose(fid)
files=cellstr(ls([dir_in, '*.tif']));
disp(files)
% fileQueue=[1:length(files)];
fileQueue=[285]; %3 for YF %285 for Sask1
exclude=[];
fileQueue=setdiff(fileQueue, exclude);
RegionGrowing=0; % set to test on global NDWI only
% tileSize has to be a multiple of 16, and apparentely
% needs to be same as processing window size

% tileSize      = [inFileInfo.TileWidth*8, inFileInfo.TileLength*8];
% tileSize      = [2048, 2048];
tileSize      = [4096, 4096];
parallel=1;
%     parpool(4);
%% Loop
for i=fileQueue
    disp(datetime)
    fprintf('File number: %d\n', i)
    name_in=files{i}; %27
    img_in=[dir_in, name_in];
    disp(dir(img_in))
    
    % Format tiff object
    inFileInfo    = imfinfo(img_in);
    
    % Parallel
%     if inFileInfo.FileSize < 3e+9
%         parallel=0;
%     else
%         parallel=0;
%     end    
    
    % Format name out
    name_out=[name_in, '_batchClass.tif'];
    img_out=[dir_out, name_out]; %NB means not border

   
    % Process images
    if RegionGrowing==1
        g = @OBIA_BP_Fun_3;
    else
        g= @OBIA_Global_BP_Fun_3;
        name_out=[name_in, '_batchClass_Global.tif'];
        img_out=[dir_out, name_out];
    end
    
    outFileWriter = BP_bigTiffWriterEK(img_out, inFileInfo(1).Height,...
        inFileInfo(1).Width, tileSize(1), tileSize(2));

    % g= @(M_strxr) M_strxr.data; % identity function as test

    window=tileSize; % block proc window size
    tic
    disp('Classifying...')
    blockproc(img_in, window, g, 'Destination', outFileWriter, 'UseParallel', parallel);
    fprintf('Done.  \n\tParallel option = %u.  Window = %u by %u pixels\n', parallel, window)
    toc

    % Save georef info
        % .mat file
    info=geotiffinfo(img_in);
%     gti_out=[dir_out, name_out(1:end-4), '.mat'];
%     save(gti_out, 'info')

        % .tfw world file
    gti_out=[dir_out, name_out(1:end-4), '.tfw'];
    worldfilewrite(info.SpatialRef, gti_out)
        % add geotiffwrite for ease (extra processing)!

    % Display
    disp('Georef files written.')
    fprintf('Output written: %s\n', img_out);
end

clear outFileWriter
outFileInfo    = imfinfo(img_out);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Batch finished.'); disp(datetime)
toc
