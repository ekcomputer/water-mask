% function to calculate area changes in lakes
clear
% close all
    % load data
pairs_path='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\pairs.mat';
load(pairs_path);
    % inuvik
% im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170716_S03X_Ch063v034_V1.tif';
% im_b_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170807_S01X_Ch063v034_V1.tif';
base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
Q=23:length(pairs); % files to load
save_geo=0; % save result?
save_mat=1; % save .mat data file
    % file formating
matPath_base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\areaChange\mat\';
matPath=[matPath_base, 'Change.mat'];
if save_mat
    load(matPath); % load summary change stats
end
% Q=30
% key: inuvik =3, yellowknife 1= 26, yellowknife 2=30 (doesn't work)

%% loop for each file pair
for n=Q % iterate over all file pairs
%     close all
    fprintf('\nFile %d/%d: %s | %s\n', pairs(n).id(1),pairs(n).id(end),...
        pairs(n).a, pairs(n).b)
    im_a_pth=[base,pairs(n).a];
    im_b_pth=[base,pairs(n).b];
    [im{1}, R_a]=geotiffread(im_a_pth);
    [im{2}, R_b]=geotiffread(im_b_pth);

    %% prep to pad smaller file to align pixels in geographic locationshelp
        % find out which file is most NW (assuming cols start from N and rows
        % start from W
    R=[R_a; R_b];
    % im={im_a, im_b}; % make into a cell array
        % find which image to pad on each side
    [top,north]=min(vertcat(R(:).YWorldLimits)); northFile=north(1); % file with most N upper extent
        top=max(top);
    [left,west]=max(vertcat(R(:).XWorldLimits)); westFile=west(1); % file with most W left extent
        left=min(left);
    [btm,south]=max(vertcat(R(:).YWorldLimits)); southFile=south(2); % file with most W left extent
        btm=min(btm);
    [right,east]=min(vertcat(R(:).XWorldLimits)); eastFile=east(2); % file with most W left extent
        right=max(right);
           
        % compute pad
    for i=1:length(im)
        [NW_intrins(i,1), NW_intrins(i,2)] = worldToIntrinsic(R(i),left, top); % outside bounds. order: x,y
        [SE_intrins(i,1), SE_intrins(i,2)] = worldToIntrinsic(R(i),right, btm);


        buf(i,1)=NW_intrins(i,2)-R(i).YIntrinsicLimits(1); % top buffer
        buf(i,2)=NW_intrins(i,1)-R(i).XIntrinsicLimits(1); % left buffer
        buf(i,3)=-SE_intrins(i,2)+R(i).YIntrinsicLimits(2); % btm buffer
        buf(i,4)=-SE_intrins(i,1)+R(i).XIntrinsicLimits(2); % right buffer

    %     buf(i,1)=-NW_intrins(i,2); % top buffer
    %     buf(i,2)=-NW_intrins(i,1); % left buffer
    %     buf(i,3)=-SE_intrins(i,2)+R(i).YIntrinsicLimits(2); % btm buffer
    %     buf(i,4)=SE_intrins(i,1)-R(i).XIntrinsicLimits(2); % right buffer

    end
    buf
    % buf=min(zeros(length(im), 4), buf)

    %% clip to align rasters % maybe use floor instead of round?
        % pad images
        buf=round(buf);
    for i=1:length(im)
%         im{i}=[-99*ones(round(buf(i,1)), size(im{i}, 2)); im{i}; -99*ones(round(buf(i,3)), size(im{i}, 2))]; % pad top and bottom
%         im{i}=[-99*ones(size(im{i}, 1), round(buf(i,2))), im{i}, -99*ones(size(im{i}, 1), round(buf(i,4)))]; % pad top and bottom
        im{i}=im{i}(1+buf(i,1):end-buf(i,3), 1+buf(i,2):end-buf(i,4));
    end


    % top-top2; add to northFile % repeat for S,E,W

    % im=im_a;
    % im(:,:,3)=im_b;

    %% stack images
    im_change=im{1}; im_change(:,:,2)= im{2};
    d.max_water=im{1}~=-99 & im{2}~=-99 & (im{1}==1 | im{2}==1); % all regionss w 2 observations that had water at least once
    % d.change=im{2}-im{1} & ~d.mask;
    d.growth =im{2}==1 & im{1}==0;
    d.shrink =im{1}==1 & im{2}==0;
%     imagesc(d.growth); axis image; title('Growth')
%     figure; imagesc(d.shrink); axis image; title('Shrink')
    

    %% stats
    disp('Calculating region properties...')
    stats=regionprops(d.max_water, 'Area', 'Perimeter', 'PixelIdxList', 'Centroid')
        % check that files actually overlap!
    %     if R(1).XWorldLImits

    if isempty(stats)
        disp('No overlap between files...skipping.')
        continue
    end
    fprintf('%d regions found.  Calculating area change...\n', length(stats))
                % create ref. object
    R_max=R(1);
        R_max.RasterSize=size(im_change(:,:,1));
        R_max.XWorldLimits=[left, right];
        R_max.YWorldLimits=[btm, top];
            % calc centroids
    centroids_intrins=vertcat(stats.Centroid);
    clear centroids
    [centroids(:,1), centroids(:,2)]=intrinsicToWorld(R_max,...
        centroids_intrins(:,1), centroids_intrins(:,2));
        % loop vars
    file_a=pairs(n).a;     
    file_b=pairs(n).b;        
    filecode_a=pairs(n).id(1);
    filecode_b=pairs(n).id(2);
    date_a=pairs(n).a(4:11);
    date_b=pairs(n).b(4:11);
    days_apart=days(datetime(date_b, 'InputFormat','yyyyMMdd')-...
        datetime(date_a, 'InputFormat','yyyyMMdd')); % number of days between observations
        % loop
    for i=1:length(stats)
        stats(i).change=sum(d.growth(stats(i).PixelIdxList))-...
            sum(d.shrink(stats(i).PixelIdxList));
        stats(i).size_a=sum(im{1}(stats(i).PixelIdxList)==1);
        stats(i).size_b=sum(im{2}(stats(i).PixelIdxList)==1);
        stats(i).pchange=stats(i).change/stats(i).size_a*100; % percent change
        gtinfo=geotiffinfo(im_a_pth);
        mstruct=geotiff2mstruct(gtinfo); % get map projection structure to convert to lat/long
        [stats(i).lat, stats(i).long] = projinv(mstruct,...
            centroids(i,1), centroids(i,2)); % add WGS84 lat/long 
        if stats(i).pchange==-100
            stats(i).pchange=-9999; % change -100 to -9999 to indicate size_b=0
        elseif isinf(stats(i).pchange)
            stats(i).pchange=9999; % in case size_a was 0
        end
        stats(i).file_a=file_a;
        stats(i).file_b=file_b;        
        stats(i).filecode_a=filecode_a;
        stats(i).filecode_b=filecode_b;
        stats(i).date_a=date_a;
        stats(i).date_b=date_b;
        stats(i).days=days_apart;
    end

    %% create shapefile before rm fields
    p= mappoint(centroids(:,1), centroids(:,2), stats);
      
        
    %% create mat file
    if save_mat
            % raster calcs
        d.mask=im{1}==-99 | im{2}==-99;

            % create new fields
        stats_all(filecode_a).change_sum.grow=sum(d.growth(:))-sum(d.shrink(:)); % save this tile's stats to all_stats file, indexed by first file's code
        stats_all(filecode_a).change_sum.valid_px=sum(d.mask(:));
            % transfer sum data
        stats_all(filecode_a).change_sum.file_a=stats(1).file_a;
        stats_all(filecode_a).change_sum.file_b=stats(1).file_b;
        stats_all(filecode_a).change_sum.filecode_a=stats(1).filecode_a;
        stats_all(filecode_a).change_sum.filecode_b=stats(1).filecode_b;
        stats_all(filecode_a).change_sum.date_a=stats(1).date_a;
        stats_all(filecode_a).change_sum.date_b=stats(1).date_b;
        stats_all(filecode_a).change_sum.days=stats(1).days;
        stats_all(filecode_a).change_sum.regions=length(stats);
        
            % remove extra data from stats
        stats=rmfield(stats, {'PixelIdxList', 'file_a', 'file_b', 'filecode_a',...
            'filecode_b', 'date_a', 'date_b', 'days'});  
            % transfer data
        stats_all(filecode_a).change=stats; % save this tile's stats to all_stats file, indexed by first file's code
            % save
        disp(['Saving .mat file to directory: ', matPath_base])
        save(matPath, 'stats_all');
    end
    if save_geo

        %% save shapefile

        shapePath_base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\areaChange\shp\';
        shapePath=[shapePath_base, 'Cnt_',pairs(n).a(18:26),'_',pairs(n).a(4:11) ,'_',pairs(n).b(4:11),  '.shp'];
        disp(['Saving shapefile to directory: ', shapePath_base])
        shapewrite(p, shapePath);

        %% save raster
        disp('Creating raster...')
        d.water=im{1}==1 & im{2}==1;
            %key: -99= no valid comparison or no data; 0=constant land; 1 =
            %constant water; 2=lost water; 3=gained water (per pixel)
        d.change_rast=int8(-99*d.mask + d.water + 2*d.shrink + 3*d.growth); % leftover values =0=constant land
        rastPath_base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\areaChange\rast\';
        rastPath=[rastPath_base, 'Chg_',pairs(n).a(18:26),'_',pairs(n).a(4:11) ,'_',pairs(n).b(4:11), '.tif'];
        disp(['Saving raster to directory: ', rastPath_base])
        geotiffwrite(rastPath, d.change_rast, R_max, 'GeoKeyDirectoryTag',gtinfo.GeoTIFFTags.GeoKeyDirectoryTag)
        imagesc(d.change_rast, [-1 3]); axis image; title('Change'); pause(0.5)
    end
    end
% end