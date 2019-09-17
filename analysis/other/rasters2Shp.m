% script to load all files and calculate size and abundance distributions

% TODO: load using matfiles? Save longitude as field
% NOTE: min size is commented out
% Note: runs fixWaterDistribuion automatically afterwards
% TODO: figure out how to subtract holes and fix topology...or remove
% references to shp1 and holes....Grrrrrrrr....
%% params
clear
addpath D:\Dropbox\Matlab\DownloadedCode\MinBoundSphere&Circle\MinBoundSphere&Circle
addpath D:\Dropbox\Matlab\DownloadedCode\MinBoundSphere&Circle\MinBoundSphere&Circle\Auxiliary
minSize=40;%40
startMid=0; % start from middle
start=1;
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
    shp_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\shp\WaterMask.shp';
    [tbl, tbl_raw]=xlsread(tbl_in);  
end


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
    msk(NoValues)=0;msk=logical(msk);
    L=bwlabel(msk);
 
        % calc regionprops
    clear init_stats
    abun(i).stats=regionprops(L>0, 'Area','Perimeter', 'Centroid','MajorAxisLength', 'PixelList');
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
        if length(abun(i).stats)==0 % if no water in image!
            [abun(i).stats.lat, abun(i).stats.long] = []; % add WGS84 lat/long 
            [abun(i).stats.SDF, abun(i).stats.LehnerDevel] = [];
            abun(i).stats(1).Area=[];
            abun(i).stats(1).Centroid=[];
            abun(i).stats(1).MajorAxisLength=[];
            abun(i).stats(1).Perimeter=[];
            abun(i).stats(1).lat=[];
            abun(i).stats(1).long=[];
            abun(i).stats(1).minBoundRadius=[];
            abun(i).stats(1).LehnerDevel=[];
            abun(i).stats(1).SDF=[];
            abun(i).freq_min40=[];
        else
        end
    
        % save as shapefile
    bound=bwboundaries(msk, 4);   
    bound_noHoles=bwboundaries(msk, 'noholes');  

%     convert=@(x) geoshape(x(:,2), x(:,1), 'Simplify', false, 'Geometry', 'polygon');
    convert=@(x) polyshape(x, 'Simplify', false);
    shp0=cellfun(convert, bound);
    shp0=simplify(shp0);
    shp1=cellfun(convert, bound_noHoles);
    shp1=simplify(shp1);
    shp_holes=xor(shp0, shp1); % <------------------here
        % retry 4
    
    shp=mapshape(); shp.Geometry='polygon';

    for j=1:length(bound)
        [lat, long]=intrinsicToWorld(R, shp0(j).Vertices(:,2), shp0(j).Vertices(:,1));
        shp(j)=mapshape(lat, long);
    end
        % delete empty features...
    j=1;
    while j <=length(shp)
        if isempty(shp(j))
            shp(j)=[];
            shp0(j)=[];
%             abun(i).stats(j)=[];
        end
        j=j+1;
    end
        shp.area = {''};
    shp.perim = {''};
    shp.centroid = {''};
    shp.sdf = {''};
%     shp.area = num2cell(area(shp0));
%     shp.perim = num2cell(perimeter(shp0));
%     shp.centroid = num2cell(centroid(shp0));
%     shp.sdf = num2cell(perimeter(shp0)./(2*sqrt(pi*area(shp0))));
%     % retry 3
%     
%     shp1=mapshape();
%     for j=1:length(bound)
%         shp1(j)=mapshape(shp(j).Vertices(:,2), shp(j).Vertices(:,1));
%     end
%     
%     % retry 2
%     
%     shp1=geoshape();
%     for j=1:length(bound)
% %         shp1(j).X=bound{j}(:,2);
% %         shp1(j).X=bound{j}(:,2);
%         [shp1(j).Latitude, shp1(j).Longitude]=intrinsicToWorld(R, shp0(j).Vertices(:,2), shp0(j).Vertices(:,1));
% %         shp1(j).Latitude(end+1)=NaN;
% %         shp1(j).Longitude(end+1)=NaN;
%     end
%     
%         % retry 1
%         
%     for j=1:length(bound)
%        [shp1(j).long, shp1(j).lat]=intrinsicToWorld(R, bound{j}(:,2), bound{j}(:,1));
%        shp(j)=geoshape(shp1(j).long, shp1(j).lat,'geometry', 'polygon');
% %        shp.SDF{j}=abun(i).stats(j).SDF;
%        [abun(j).longitude, abun(j).latitude]=intrinsicToWorld(R, bound{j}(:,2), bound{j}(:,1));
% %        shp(j)=mapshape(abun(i).stats(j), 'Geometry', 'polygon');
%     end
%     shp.Geometry='polygon';
%     foo=geoshape(rmfield(abun(i).stats, {'lat','long'}))
%     %, 'Geometry', 'polygon');
    
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
        shapewrite(shp, shp_out)
        fprintf('Saving to %s\n', shp_out)
    end
end
shapewrite(shp, shp_out)
fprintf('Saving to %s\n', shp_out)
disp('Done.')
