%BP_OBIA_devel uses new version of classification written for block
%processing, but allows user to implement on small regions w/o block
%processing (for development purposes)

% todo: add geotiffwrite at end insted of .tfw; change bigtiff to tiff
% writer; edit optomizeConn_16 for case of no peak in MasterMetric and
% other scaling issues (incl no bimodal hist...)
% Script to apply block processing to DCS images- for water classification
% modified from DCS_images_4.m

% File queue
clear; close all; clc
tic
% set(0,'DefaultFigureVisible','off')
% dir_in='D:\ArcGIS\FromMatlab\ClipSquares\';
% dir_in='D:\GoogleDrive\ABoVE top level folder\AirSWOT_CIR\DAAC_Preview\DCS_all\';
dir_in='D:\ArcGIS\FromMatlab\ClipSquares\';
dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\';


f.logDir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\logs\';
logfile=[f.logDir, 'log.txt'];


% logfile='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\logs\log.txt';
% fid=fopen(logfile, 'a');
% fprintf(fid, '----------------------\n');
% fclose(fid);
files=cellstr(ls([dir_in, '*.tif']));
disp('Files:')
disp([num2cell([1:length(files)]'), files])

% fileQueue=[1:length(files)];
fileQueue=[16]; %3 for YF %285 for Sask1
% fileQueue=[29:-1:2]; %3 for YF %285 for Sask1

% fileQueue=find(files=="LakeChange0807.tif");
% fileQueue=find(files=="cs_14_RedberCorner.tif");
exclude=[];
fileQueue=setdiff(fileQueue, exclude);
RegionGrowing=1; % set to test on global NDWI only
% tileSize has to be a multiple of 16, and apparentely
% needs to be same as processing window size

parallel=0;
%     parpool(4);
datecode=char(datetime('now','Format','yyyy-MM-dd-HHmm'));

disp(datetime)
%% Loop
for i=fileQueue
    fprintf('File number: %d\n', i)
    name_in=files{i}; %27
    img_in=[dir_in, name_in];
    fprintf('Classifying file:\t%s\n', name_in);
    [cir, R]=geotiffread(img_in);

    % Process images
    tic
    disp('Classifying...')
    if RegionGrowing==1
        name_out=[name_in(1:end-4), '_batchClass.tif'];
        classified_out=OBIA_BP_Fun(cir, f.logDir, 'local', name_out, datecode);
        img_out=[dir_out, name_out]; %NB means not border
    else
        name_out=[name_in(1:end-4), '_batchClass_Global.tif'];
        classified_out=OBIA_BP_Fun(cir, f.logDir, 'global', name_out, datecode);
        img_out=[dir_out, name_out];
    end
    fprintf('Done.  Run in developer (no block proc) mode.\n')

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

    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    disp('Class test finished.'); disp(datetime)
    elapsedTime=toc; fprintf('Elapsed time:\t%3.2f minutes\n', toc/60);

    %% Save
    try
    fprintf('Saving to directory: %s\n', dir_out)
    geotiffwrite(img_out, classified_out,...
        R, 'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
    catch
        disp('File not saved.')
    end
    
    %% for saving plots
    
    figure(6)
end