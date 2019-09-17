for i=1:length(files)
    fprintf('%s\n', files{i})
    paths{i}=[dir_in, files{i}];
    paths_out{i}=[dir_out, files_out{i}];
end

%%

files_out=cellstr(ls([dir_out, '*.tif']))
for i=1:length(files_out)
    fprintf('%s\n', files_out{i})
end
%% Read table

tbl_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\WC_LOG.xlsx';
[tbl, tbl_raw]=xlsread(tbl_in, 1);

%% practice visulizing
for i=6:length(paths_out)
    fprintf('File %d:\t%s', i, files_out{i})
    WC=imread(paths_out{i});
%     imshow(paths_out{i}, 'DisplayRange', [0 1]); %colormap parula
    WC(WC==-99)=-1;
    figure(1); imagesc(WC); axis image
%     figure(2); imshow(paths{i})
    fprintf('\n')
    pause
end