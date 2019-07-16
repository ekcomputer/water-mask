% script to move all files with unique AS coverage
% remember to manually move straggler

%% params
clear
minSize=40;%40
startMid=0; % start from middle
start=1;
atUCLA=0;
copyShape=1;
copyTiff=0;
%% directories
if atUCLA
    mask_dir='J:\Final\';
    file_list_in='J:\Final\logs\WC_LOG_Summ.xlsx';
    tbl_in='J:\output\analysis\total_list2.xlsx';
    struct_out='J:\output\analysis\distrib.mat';
    [tbl, tbl_raw]=xlsread(tbl_in);
else % at brown
    mask_dir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
    file_list_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\WC_LOG_Summ.xlsx';
    tbl_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\total_list2.xlsx';
    struct_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\distrib.mat';
    [tbl, tbl_raw]=xlsread(tbl_in);  
    in_shp_dir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_simpl\';
    out_shp_dir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_simpl_unq\';
    out_tiff_dir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final_unq\';
end
% [~, files]=xlsread(file_list_in, 1);
% list=files(:,1);

%% load input data
all_files=cellstr(ls([mask_dir, '*.tif']));
use=tbl(:,4); use=use==1 | use==2;
if startMid
    load(struct_out);
end
files=all_files(use);
%% loop

for i= start:length(files) %+1
    disp(i)
    if i==length(files)+1
       disp('Adding straggler.')
       pth_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\straggler\WC_20170806_S01X_Ch045v030_V1_clip.tif';
    end
    if copyShape
        in_shp={[in_shp_dir, files{i}(1:end-4), '.shp'],
            [in_shp_dir, files{i}(1:end-4), '.dbf'],
            [in_shp_dir, files{i}(1:end-4), '.prj'],
            [in_shp_dir, files{i}(1:end-4), '.shx']};
        out_shp={[out_shp_dir, files{i}(1:end-4), '.shp'],
            [out_shp_dir, files{i}(1:end-4), '.dbf'],
            [out_shp_dir, files{i}(1:end-4), '.prj'],
            [out_shp_dir, files{i}(1:end-4), '.shx']};
        for j=1:4
            copyfile(in_shp{j}, out_shp{j})
        end
    end
    if copyTiff
        in_tiff=[mask_dir, files{i}];
        out_tiff=[out_tiff_dir, files{i}];
        copyfile(in_tiff, out_tiff)
    end
end