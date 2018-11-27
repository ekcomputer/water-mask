% function to calculate area changes in lakes
clear
    % load data
    % inuvik
im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170716_S03X_Ch063v034_V1.tif';
im_b_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170807_S01X_Ch063v034_V1.tif';
%     % yellowknife
% im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170715_S02X_Ch083v065_V1.tif';
% im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170812_S03X_Ch083v065_V1.tif';
%     % yellowknife 2
% im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170715_S03B_Ch084v067_V1.tif';
% im_a_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\WC_20170812_S04X_Ch084v067_V1.tif';

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
    buf(i,3)=-SE_intrins(i,2)+R(i).YIntrinsicLimits(2); % btm buffer
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

%%
im_change=im{1}; im_change(:,:,2)= im{2};
% d.mask=im{1}==-99 | im{2}==-99;
d.max_water=im{1}==1 & im{2}==1;
% d.change=im{2}-im{1} & ~d.mask;
d.growth =im{2}==1 & im{1}==0;
d.shrink =im{1}==1 & im{2}==0;

stats=regionprops(d.max_water, 'Area', 'Perimeter', 'PixelIdxList', 'Centroid')