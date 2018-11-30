% function to calculate area changes in lakes
clear
close all
    % load data
pairs_path='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\pairs.mat';
load(pairs_path);
    % inuvik
% im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170716_S03X_Ch063v034_V1.tif';
% im_b_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170807_S01X_Ch063v034_V1.tif';
base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
Q=1:length(pairs); % files to load
save=1; % save result?
Q=[3]
% key: inuvik =3, yellowknife 1= 26, yellowknife 2=30 (doesn't work)

%% loop for each file pair
for n=Q % iterate over all file pairs
    fprintf('File %d/%d: %s | %s\n', pairs(n).id(1),pairs(n).id(end),...
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
    [top,north]=max(vertcat(R(:).YWorldLimits)); northFile=north(1); % file with most N upper extent
        top=max(top);
    [left,west]=min(vertcat(R(:).XWorldLimits)); westFile=west(1); % file with most W left extent
        left=min(left);
    [btm,south]=min(vertcat(R(:).YWorldLimits)); southFile=south(2); % file with most W left extent
        btm=min(btm);
    [right,east]=max(vertcat(R(:).XWorldLimits)); eastFile=east(2); % file with most W left extent
        right=max(right);
           
        % compute pad
    for i=1:length(im)
        [NW_intrins(i,1), NW_intrins(i,2)] = worldToIntrinsic(R(i),left, top); % outside bounds. order: x,y
        [SE_intrins(i,1), SE_intrins(i,2)] = worldToIntrinsic(R(i),right, btm);


        buf(i,1)=-NW_intrins(i,2)+R(i).YIntrinsicLimits(1); % top buffer
        buf(i,2)=-NW_intrins(i,1)+R(i).XIntrinsicLimits(1); % left buffer
        buf(i,3)=SE_intrins(i,2)-R(i).YIntrinsicLimits(2); % btm buffer
        buf(i,4)=SE_intrins(i,1)-R(i).XIntrinsicLimits(2); % right buffer

    %     buf(i,1)=-NW_intrins(i,2); % top buffer
    %     buf(i,2)=-NW_intrins(i,1); % left buffer
    %     buf(i,3)=-SE_intrins(i,2)+R(i).YIntrinsicLimits(2); % btm buffer
    %     buf(i,4)=SE_intrins(i,1)-R(i).XIntrinsicLimits(2); % right buffer

    end
    buf
    % buf=min(zeros(length(im), 4), buf)

    %% pad % maybe use floor instead of round?
        % pad images
    for i=1:length(im)
        im{i}=[-99*ones(round(buf(i,1)), size(im{i}, 2)); im{i}; -99*ones(round(buf(i,3)), size(im{i}, 2))]; % pad top and bottom
        im{i}=[-99*ones(size(im{i}, 1), round(buf(i,2))), im{i}, -99*ones(size(im{i}, 1), round(buf(i,4)))]; % pad top and bottom
    end


    % top-top2; add to northFile % repeat for S,E,W

    % im=im_a;
    % im(:,:,3)=im_b;

    %% stack images
    im_change=im{1}; im_change(:,:,2)= im{2};
    % d.mask=im{1}==-99 | im{2}==-99;
    d.max_water=im{1}~=-99 & im{2}~=-99 & (im{1}==1 | im{2}==1); % all regionss w 2 observations that had water at least once
    % d.change=im{2}-im{1} & ~d.mask;
    d.growth =im{2}==1 & im{1}==0;
    d.shrink =im{1}==1 & im{2}==0;
    imagesc(d.growth); axis image; title('Growth')
    figure; imagesc(d.shrink); axis image; title('Shrink')

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
    for i=1:length(stats)
%         disp(i)
        stats(i).change=sum(d.growth(stats(i).PixelIdxList))-...
            sum(d.shrink(stats(i).PixelIdxList));
        stats(i).size_a=sum(im{1}(stats(i).PixelIdxList)==1);
        stats(i).size_b=sum(im{2}(stats(i).PixelIdxList)==1);
        stats(i).pchange=stats(i).change/stats(i).size_a*100; % percent change
        if stats(i).pchange==-100
            stats(i).pchange=-9999; % change -100 to -9999 to indicate size_b=0
        elseif isinf(stats(i).pchange)
            stats(i).pchange=9999; % in case size_a was 0
        end
    end

    %% create shapefile
        % create ref. object
    R_max=R(1);
        R_max.RasterSize=size(im_change(:,:,1));
        R_max.XWorldLimits=[left, right];
        R_max.YWorldLimits=[btm, top];

        %format shapefile
    centroids_intrins=vertcat(stats.Centroid);
    [centroids(:,1), centroids(:,2)]=intrinsicToWorld(R_max, centroids_intrins(:,1), centroids_intrins(:,2));
    p= mappoint(centroids(:,1), centroids(:,2), stats);

    if save
        %% save shapefile

        shapePath_base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\areaChange\shp\';
        shapePath=[shapePath_base, 'Cnt_',pairs(n).a(18:26),'_',pairs(n).a(4:11) ,'_',pairs(n).b(4:11),  '.shp'];
        disp(['Saving shapefile to directory: ', shapePath_base])
        shapewrite(p, shapePath);

        %% save raster
        rastPath_base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\areaChange\rast\';
        rastPath=[rastPath_base, 'Chg_',pairs(n).a(18:26),'_',pairs(n).a(4:11) ,'_',pairs(n).b(4:11), '.tif'];
            % this is slow.  use matrix element mult.?
        d.change_rast=2*(d.growth-d.shrink); % growth is 1, shrink is -1
        d.change_rast(im{1}==-99 | im{2}==-99)=-99; % mark NoValue areas
        d.change_rast(im{1}==1 & im{2}==1)=1; % mark constant water
        d.change_rast(im{1}==0 & im{2}==0)=0; % mark constant land       
        gtinfo=geotiffinfo(im_a_pth);
        disp(['Saving raster to directory: ', rastPath_base])
        geotiffwrite(rastPath, d.change_rast, R_max, 'GeoKeyDirectoryTag',gtinfo.GeoTIFFTags.GeoKeyDirectoryTag)
    end
end