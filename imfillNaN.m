function filled=imfillNaN(bw, mask)

% V2 is faster
% filled=imfillNaN(bw, mask)
% Function to fill in NaN values surrounded by foreground
% in classified, binary image
% mask is binary mask with 1= no value (NaN value)


%% 
origPixels=sum(sum(bw));
% cc1=bwconncomp(~bw); cc1.PixelIdxList=[];
holes = imfill(bw,'holes') &~bw; % can include real islands
holes=holes&mask; % re using name to save space.  holes now excludes real isalnds
filled=bw;
filled(holes==1)=1;
finalPix=sum(sum(filled));
% cc2=bwconncomp(~filled); cc2.PixelIdxList=[];
cc=regionprops(holes, 'Area');
fprintf('\tFilled all NoData cells completely surrounded by water class.\n')
fprintf('\tNumber of cells filled:\t%d\n', finalPix-origPixels)
fprintf('\tNumber of regions filled:\t%i\n', length(cc))
%% fill



% filled=logical(round(regionfill(uint8(bw), ~mask)));