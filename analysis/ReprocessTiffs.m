% for re-writing tiffs from corrupted mask (WM) files!
% uses georeff info from non-corrupt image files.
%%i/o
clear
% dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
dir_final='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';

files_in=cellstr(ls([dir_final, '*.tif'])); % this is wrong directory, but it will give me correct filenames for IDS in excel sheet

% dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final2\';
dir_out=dir_in;

dir_georef='F:\AboveDCSRasterManagement\CanadaAlbersTranslate\';
georef_files=cellstr(ls([dir_georef, '*.tif']));

    % make queue
% Q=17:29;
    % 11/15 REPROCESS These are file number from 'final' folder
% Q=[35	36	37	38	39	40	42	43	44	56	57	58	61];
    % 11/16 REPROCESS These are file number from 'final' folder
% Q=[44	62	63	64	66	67	69	71	74	75	76	77	78	79	80	82	85	87	88	89	90	91	92	93	94	95	96	97	98	99	105	106	107	108	112	113	114	116	121	124	127	129	130	132	136	142	146	147	148	151	152	154	155	161	166	167	179	180	181	182	185	186	187	191	193	198	201	203	204	206	207	210	211	219	];
    % 11/18: just need to run 2
Q=[96, 98];
disp('Files in Q:')
(files_in{Q})
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