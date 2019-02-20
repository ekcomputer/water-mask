% script to load all files and calculate size and abundance distributions

% TODO: load using matfiles? Save longitude as field
% NOTE: min size is commented out
% Note: runs fixWaterDistribuion automatically afterwards
%% params
clear
addpath D:\Dropbox\Matlab\DownloadedCode\MinBoundSphere&Circle\MinBoundSphere&Circle
addpath D:\Dropbox\Matlab\DownloadedCode\MinBoundSphere&Circle\MinBoundSphere&Circle\Auxiliary
minSize=40;%40
startMid=1; % start from middle
start=251;
atUCLA=0;
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
end
% [~, files]=xlsread(file_list_in, 1);
% list=files(:,1);

%% load input data
all_files=cellstr(ls([mask_dir, '*.tif']));
use=tbl(:,4); use=use==1 | use==2;
if startMid
    load(struct_out);
end
%% loop
files=all_files(use);
for i= start:length(files)+1
    disp(i)
    if i==length(files)+1
       disp('Adding straggler.')
       pth_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\straggler\WC_20170806_S01X_Ch045v030_V1_clip.tif';
    else
        pth_in=[mask_dir, files{i}];
    end
    [msk, R]=geotiffread(pth_in);
    NoValues=(msk==-99);
        % calc limnicity (water fraction)
    abun(i).lim=sum(sum(msk==1))/sum(sum(msk>=0));
    abun(i).land=sum(sum(msk==0));
    abun(i).water=sum(sum(msk==1));  
        % buffer all water region
    msk(NoValues)=0;msk=logical(msk);
    L=bwlabel(msk);
    init_stats=regionprops(L, 'PixelIdxList', 'Area');
    SE=strel('diamond',1);
    NoValues_dil=imdilate(NoValues, SE); % msk is now dilated by one
    NoValues_dil(:,1)=1; NoValues_dil(:,end)=1;
    NoValues_dil(1,:)=1; NoValues_dil(end,:)=1;
    
        % remove edge water
    edgeRegions=setdiff(unique(L(NoValues_dil)), 0);
    edgeRegionPx=vertcat(init_stats(edgeRegions).PixelIdxList);
    L(edgeRegionPx)=0;
    
        % size filter
%     smallRegions=[init_stats.Area]<minSize;
%     smallRegionsPx=vertcat(init_stats(smallRegions).PixelIdxList);
%     L(smallRegionsPx)=0;
 
        % calc regionprops
    clear init_stats
    abun(i).stats=regionprops(L>0, 'Area','Perimeter', 'Centroid','MajorAxisLength', 'PixelList');
    if length(abun(i).stats)==0 % if no water in image!
        [abun(i).stats.LehnerDevel, abun(i).stats.long] = []; % add WGS84 lat/long 
        [abun(i).stats.SDF, abun(i).stats.LehnerDevel] = [];
    else
            % plot histogram of stats
        h=histogram([abun(i).stats.Area]/1e6); xlabel('Area ($km^2$)'); ylabel('Count');
        if i<=length(files)
        title(files{i}, 'Interpreter', 'none')
        end
    %         set(gca, 'XScale', 'log', 'YScale', 'log')
        drawnow
        intrinsic=vertcat(abun(i).stats.Centroid);
        clear world
        [world(:,1), world(:,2)]=intrinsicToWorld(R,intrinsic(:,1), intrinsic(:,2));
        gtinfo=geotiffinfo(pth_in);
        mstruct=geotiff2mstruct(gtinfo); % get map projection structure to convert to lat/long
        for j=1:length(abun(i).stats)
            [abun(i).stats(j).lat, abun(i).stats(j).long] = projinv(mstruct,...
                world(j,1), world(j,2)); % add WGS84 lat/long 
            abun(i).stats(j).minBoundRadius=ExactMinBoundCircle(abun(i).stats(j).PixelList);
            abun(i).stats(j).LehnerDevel=(abun(i).stats(j).Area)/(abun(i).stats(j).minBoundRadius^2*pi);
            abun(i).stats(j).SDF=abun(i).stats(j).Perimeter/(2*sqrt(pi*abun(i).stats(j).Area));
        end
        abun(i).freq_min40=length(abun(i).stats)/(abun(i).water+abun(i).land)*100e6;
        abun(i).stats=rmfield(abun(i).stats, 'PixelList');
    end
    
        % save as shapefile
%     bound=bwboundaries(msk);   
%     bound_noHoles=bwboundaries(msk, 'noholes');   
%     convert=@(x) polyshape(x, 'Simplify', true);
%     shp=cellfun(convert, bound);
        % add addiitonal data
    if i~=length(files)+1
        abun(i).file=files{i};
        abun(i).file_idx=find(strcmp(all_files, abun(i).file));
    else
        abun(i).file='WC_20170806_S01X_Ch045v030_V1_clip.tif';
        abun(i).file_idx=-99;
    end
        % save every ten times
    if mod(i, 10)==0
        save(struct_out, 'abun')
        fprintf('Saving to %s\n', struct_out)
    end
end
save(struct_out, 'abun')
fprintf('Saving to %s\n', struct_out)
disp('Done.')

%% Fix entries with no water
fixWaterDistribution