% todo: add geotiffwrite at end insted of .tfw; change bigtiff to tiff
% writer; edit optomizeConn_16 for case of no peak in MasterMetric and
% other scaling issues (incl no bimodal hist...)
% Script to apply block processing to DCS images- for water classification
% modified from DCS_images_4.m

% File queue
clear; clc; close all
disp('Batch started.'); disp(datetime)
tbatch=tic;
global f
dbstop if error
opengl software % don't use AMD graphics card/driver
% f.ETHANTEST='yeah!';
f.plot=false;
f.logDir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\';
% set(0,'DefaultFigureVisible','off')
dir_in='F:\AboveDCSRasterManagement\CanadaAlbersTranslate\';
dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
logfile=[f.logDir, 'log.txt'];
% fid=fopen(logfile, 'a');
% fprintf(fid, '----------------------\n');
% fclose(fid);
% files=cellstr(ls([dir_in, '*.tif']));
files=cellstr(ls([dir_in, '*.tif']));
disp('Files:')
disp([num2cell([1:length(files)]'), files])
% fileQueue=[1:length(files)];
% fileQueue=[285]; %3 for YF %285 for Sask1 %30 for PAD
exclude=[];
% fileQueue=find(files=="DCS_20170709_S02X_Ch081v092_V1.tif");
% fileQueue=find(files=="DCS_20170716_S01X_Ch066v032_V1.tif");
% fileQueue=find(files=="DCS_20170716_S01X_Ch066v033_V1.tif");
% fileQueue=find(files=="DCS_20170709_S03B_Ch081v102_V1.tif");
% fileQueue=find(files=="DCS_20170716_S02X_Ch066v032_V1.tif"); % TK lakes w clouds
    %Testing Tiles File Queue
% fileQueue=1+[17	30	54	69	94	111	123	144	159	191	203	221	245	279	286	262 304	309	322];
% fileQueue=[18];
    % first rerun 11/13/2018: Going to Final3
% fileQueue=[3	4	11	12	14	15	16	17	18	19	23	24	26	27	28	30	31	32	33	34	35	36	37	38	39	40	42	43	44	56];
    % second rerun 11/14/2018: includes first and second time reruns: going
    % to  final 3
% fileQueue=[3	11	12	15	16	18	19	23	24	26	27	31	32	33	34	35	36	37	38	39	40	42	43	44	56	57	58	61	62	63	64	66	67	69	71	74	75	76	77	78	79	80	82	85	87	88	89	90	91	92	93	94	95	96	97	98	99	105]

    % third rerun 11/15/2018: going to final3 (overwrite)
% fileQueue=[3	11	15	16	23	24	26	27	32	34	36	37	38	39	42	43	44	62	63	64	66	67	69	71	74	75	76	77	78	79	80	82	85	87	88	89	90	91	92	93	94	95	96	97	98	99	105	106	107	108	112	113	114	116	121	124	127	129	130	132	136	142	146	147	148	151	152	154	155	161	166	167	179	180	181	182	185	186	187	191	193	198	201	203	204	206	207	210	211	219];
    % fourth rerun 11/16/2018: going to *final* (overwrite)
% fileQueue=[11	16	32	34	36	43	66	67	69	71	79	80	82	88	90	91	96	98	99	105	106	113	114	127	130	132	136	146	147	154	166	179	180	181	206	207	210	222	223	224	225	229	232	233	236	237	241	242	255	257	258	260	261];
    % fifth rerun 11/16: running new files first to prevent/delay graphics
    % erros bc small window size:
% fileQueue=[99	105	106	113	114	127	130	132	136	146	147	154	166	179	180	181	206	207	210	222	223	224	225	229	232	233	236	237	241	242	255	257	258	260	261	271	282	286	288	289	293	294	295	296	299	300	302	304	305	306	307	310	311	312	315	316	317	319	320	321	322	323	324	326	327	329	330,... % new
%     69	79	80	82	90	91	96	98,... % old
%     11	16	32	34	36	297	308	309	318]; % IR
    %11/20
% fileQueue=[11]; % test run
    %11/20 during day after graphics fail
% fileQueue=[11 16 32	34	36	69	79	80	82	90	91	98	99	106	127	146	166	179	207	225	229	233	236	237	258	271	286	289	295	297	300	302	304	306	308	309	310	311	312	316	318	319	320	321	322	323	326	327	329	330];
% fileQueue=[69	79	80	82	90	91	98	99	106	127	146	166	179	207	225	229	233	236	237	258	271	286	289	295	297	300	302	304	306	308	309	310	311	312	316	318	319	320	321	322	323	326	327	329	330];
    % 11/21 after running 3/4 during the day:
% fileQueue=[32	34	79	98	106	146	166	179	207	229	237	271	289	295	297	300	302	304	306	308	309	310	311	312	316	318	319	320	321	322	323	326	327	329	330];
    % 11/25 partial run
% fileQueue=[79	146	166	179	207	229	237	271	289	295	300];
    % 11/26
% fileQueue=[79	179	207	229	237	271	295	300	308	310	311	318	319	320	321	322	326	329	330];
    % 11/30
% fileQueue=[79	179	207	229	237	271	300	308	310	311	318	321	326	330];
    % 12/4
% fileQueue=[331 310 311 321 326 329];
    % 12/19
% fileQueue=[2	3	9	16	27	31	33	38	42	116	117	120	121	195	302];
    % 12/20
fileQueue=[116 120];
fileQueue=setdiff(fileQueue, exclude, 'stable');

RegionGrowing=1; % set to test on global NDWI only
% tileSize has to be a multiple of 16, and apparentely
% needs to be same as processing window size

% tileSize      = [inFileInfo.TileWidth*8, inFileInfo.TileLength*8];
% tileSize      = [2048, 2048];
% tileSize      = [4096, 4096];
% tileSize      = [5760, 5760];
% tileSize      = [8192, 8192];
% tileSize      = [8192, 4096];
% tileSize      = [8192, 8192]; %[5760, 5760];
% tileSize      = [16384 4096]; %4:1 aspect
% tileSize      = [12560 5376]; %2.34:1 aspect
% tileSize      = [13200 6016]; %2.19:1 aspect
% tileSize      = [16000 8000]; %2:1 aspect a little bit TOO BIG prod = 128M
    % use
% tileSize      = [9600 9600]; %1:1 aspect GOOD prod =92M
% tileSize      = [14720, 7360]; %2:1 aspect with prod 108M
    % after graphics card issues:
    % try tilesize  7200       14400 (2:1 aspect, prod=104M)
    %       6400       12800 (2:1, prod= 83M)
    %       9120        9120 (1:1, prod = 83M)
    


parallel=0;
if parallel==1
%     tileSize      = [6400, 6400];
    tileSize      = [5760, 5760];
    
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
    name_out=['WC', name_in(4:end-4), '.tif'];
    img_out=[dir_out, name_out]; %NB means not border
   
    % load params from check file
    try
        f.aConn=tbl(i,26); % 
        f.bConn=tbl(i,27); %
        tileSize=[tbl(i,28), tbl(i,29)];
        f.bounds=[tbl(i,30), tbl(i,31)];
        f.NDWIWaterAmount=tbl(i,33); %                                                                                   -
        f.NDWILandAmount=tbl(i,32);
        f.wp=tbl(i,34);
        f.windex=tbl_raw{i+1, 35}; % note raw table input includes column headers...
        f.Tlim=tbl(i,39);
    catch % just in case parsing problem
        warning('Trouble reading f params.')
        f.aConn=15; % 
        f.bConn=170; %
        tileSize=[9600 9600];
        f.bounds=[1.5 2.5];
        f.NDWIWaterAmount=0.04; %                                                                                   -
        f.NDWILandAmount=-0.06;
        f.wp=10;
        f.windex='NDWI';
        f.Tlim=5.3;
    end
    % Process images
    if RegionGrowing==1
        g = @(block_struct)  OBIA_BP_Fun(block_struct, f.logDir, 'local', img_out, datecode);
    else
        g = @(block_struct)  OBIA_BP_Fun(block_struct, f.logDir, 'global', img_out, datecode);
        name_out=['WC', name_in(4:end-4), '_Global.tif'];
        img_out=[dir_out, name_out];
    end
    
    outFileWriter = BP_bigTiffWriterEK(img_out, inFileInfo(1).Height,...
        inFileInfo(1).Width, tileSize(1), tileSize(2));

    % g= @(M_strxr) M_strxr.data; % identity function as test

    window=tileSize; % block proc window size
    tic
    disp('Classifying...')
    try
        blockproc(img_in, window, g, 'Destination', outFileWriter, 'UseParallel', parallel);
    catch
        warning('Blockproc failed.  Trying again with smaller tile size.')
        tileSize      = [8192, 8192];
        outFileWriter = BP_bigTiffWriterEK(img_out, inFileInfo(1).Height,...
            inFileInfo(1).Width, tileSize(1), tileSize(2));
        window=tileSize; % block proc window size
        blockproc(img_in, window, g, 'Destination', outFileWriter, 'UseParallel', parallel);
    end
    fprintf('Done.  \n\tParallel option = %u.  Window = %u by %u pixels\n', parallel, window)
    toc

    % Save georef info
        % .mat file
    try
        info=geotiffinfo(img_in);
    %     gti_out=[dir_out, name_out(1:end-4), '.mat'];
    %     save(gti_out, 'info')

            % .tfw world file
        gti_out=[dir_out, name_out(1:end-4), '.tfw'];
        worldfilewrite(info.SpatialRef, gti_out)
            % add geotiffwrite for ease (extra processing)!

            % .proj file
         AlbersProj='prj\Canada Albers Equal Area Conic.prj';
         proj_out=[dir_out, name_out(1:end-4), '.prj'];
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
