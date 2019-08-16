% Script to apply block processing to water classification.  Had some bugs
% with Matlab r2018a and AMD Radeon Pro Graphics card driver- causing a
% memory error and the blue screen of death.  Bug is now possibly fixed by 
% updating AMD driver.  Workaround at the time was to use opengl for
% graphics operations rather than AMD driver.  BP stands for block
% processing and OBIA stands for object-based image classification.

% File queue
clear; clc; close all
disp('Batch started.'); disp(datetime)
tbatch=tic;
global f % this structure holds all user-dependent parameters, such as file paths and tuning parameters
dbstop if error
opengl software % don't use AMD graphics card/driver
%%%%%%%%%%%%%%%%%%%% user params
exclude=[]; % exclude files from queue
fileQueue=[116 120]; % list of file numbers in input directory to process.
f.plot=false; %set as true if you want to asee plots (slows it down a little)
f.logDir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\';
% set(0,'DefaultFigureVisible','off')
f.dir_in='F:\AboveDCSRasterManagement\CanadaAlbersTranslate\';
f.dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
f.parallel=0;
f.RegionGrowing=1; % set to test on global NDWI only

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%logfile=[f.logDir, 'log.txt'];
% fid=fopen(logfile, 'a');
% fprintf(fid, '----------------------\n');
% fclose(fid);
% files=cellstr(ls([f.dir_in, '*.tif']));
files=cellstr(ls([f.dir_in, '*.tif']));
disp('Files:')
disp([num2cell([1:length(files)]'), files])
fileQueue=setdiff(fileQueue, exclude, 'stable');

% f.tileSize has to be a multiple of 16, and apparentely
% needs to be same as processing window size
 
if f.parallel==1
%     f.tileSize      = [6400, 6400];
    f.tileSize      = [5760, 5760];
    
    try parpool(2);
    catch
    end
end

datecode=char(datetime('now','Format','yyyy-MM-dd-HHmm'));
disp('File queue:')
disp(files(fileQueue))
    % read param table
tbl_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\WC_LOG_Summ.xlsx';
[tbl, tbl_raw]=xlsread(tbl_in, 1);
%% Loop
for i=fileQueue
    disp(datetime)
    fprintf('File number: %d\n', i)
    name_in=files{i}; %27
    img_in=[f.dir_in, name_in];
    disp(dir(img_in))
    
    % Format tiff object
    inFileInfo    = imfinfo(img_in);
    
    % Parallel
%     if inFileInfo.FileSize < 3e+9
%         f.parallel=0;
%     else
%         f.parallel=0;
%     end    
    
    % Format name out
    name_out=['WC', name_in(4:end-4), '.tif'];
    img_out=[f.dir_out, name_out]; %NB means not border
   
    % load params from check file
    try
        f.aConn=tbl(i,26); % 
        f.bConn=tbl(i,27); %
        f.tileSize=[tbl(i,28), tbl(i,29)];
        f.bounds=[tbl(i,30), tbl(i,31)];
        f.NDWIWaterAmount=tbl(i,33); %                                                                                   -
        f.NDWILandAmount=tbl(i,32);
        f.wp=tbl(i,34);
        f.windex=tbl_raw{i+1, 35}; % note raw table input includes column headers...
        f.Tlim=tbl(i,39);
    catch % just in case parsing problem
        warning('Trouble reading user params (structure f).')
        f.aConn=15; % 
        f.bConn=170; %
        f.tileSize=[9600 9600];
        f.bounds=[1.5 2.5];
        f.NDWIWaterAmount=0.04; %                                                                                   -
        f.NDWILandAmount=-0.06;
        f.wp=10;
        f.windex='NDWI';
        f.Tlim=5.3;
    end
    % Process images
    if f.RegionGrowing==1
        g = @(block_struct)  OBIA_BP_Fun(block_struct, f.logDir, 'local', img_out, datecode);
    else
        g = @(block_struct)  OBIA_BP_Fun(block_struct, f.logDir, 'global', img_out, datecode);
        name_out=['WC', name_in(4:end-4), '_Global.tif'];
        img_out=[f.dir_out, name_out];
    end
    
    outFileWriter = BP_bigTiffWriterEK(img_out, inFileInfo(1).Height,...
        inFileInfo(1).Width, f.tileSize(1), f.tileSize(2));

    % g= @(M_strxr) M_strxr.data; % identity function as test

    window=f.tileSize; % block proc window size
    tic
    disp('Classifying...')
    try
        blockproc(img_in, window, g, 'Destination', outFileWriter, 'UseParallel', f.parallel);
    catch
        warning('Blockproc failed.  Trying again with smaller tile size.')
        f.tileSize      = [8192, 8192];
        outFileWriter = BP_bigTiffWriterEK(img_out, inFileInfo(1).Height,...
            inFileInfo(1).Width, f.tileSize(1), f.tileSize(2));
        window=f.tileSize; % block proc window size
        blockproc(img_in, window, g, 'Destination', outFileWriter, 'UseParallel', f.parallel);
    end
    fprintf('Done.  \n\tParallel option = %u.  Window = %u by %u pixels\n', f.parallel, window)
    toc

    % Save georef info
        % .mat file
    try
        info=geotiffinfo(img_in);
    %     gti_out=[f.dir_out, name_out(1:end-4), '.mat'];
    %     save(gti_out, 'info')

            % .tfw world file
        gti_out=[f.dir_out, name_out(1:end-4), '.tfw'];
        worldfilewrite(info.SpatialRef, gti_out)
            % add geotiffwrite for ease (extra processing)!

            % .proj file
         AlbersProj='prj\Canada Albers Equal Area Conic.prj';
         proj_out=[f.dir_out, name_out(1:end-4), '.prj'];
%          copyfile(AlbersProj, proj_out);
            
            
        % Display
        disp('Georef files written.')
    catch
        disp('NO georef files written.')
    end
    fprintf('Output written: %s\n', img_out);
end

clear outFileWriter
outFileInfo    = imfinfo(img_out);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('Batch finished.'); disp(datetime)
elapsedTime=toc(tbatch); fprintf('Elapsed time:\t%3.2f minutes\n', ...
    elapsedTime/60);
