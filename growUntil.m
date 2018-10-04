function filled=imfillNaN(bw, mask)

% Function to fill in NaN values in classified, binary image
% Only fills holes that are completely surrounded by same class (or with
% max of two non-class pixels)
% mask is binary mask with 1= no value (NaN value)

%%%%%%%%%%%%%%%%%%%%%%%%testing
% clear
% close all
% dir_in='D:\ArcGIS\FromMatlab\ClipSquares\';
% disp('Clips in directory:')
% tifs=cellstr(ls([dir_in, '*.tif']))
% name_in=tifs{1}
% img_in=[dir_in, name_in]
% [cir, R]=geotiffread(img_in);
% 
% %% Set NoValue cells ==0
% 
% NoValues=or(cir(:,:,1)==65535 & cir(:,:,2)==65535 & cir(:,:,3)==65535,...
%     cir(:,:,1)==0 & cir(:,:,2)==0 & cir(:,:,3)==0); %255
% if strcmp(class(cir),'single')
%    cir=cir.*repmat(single(~NoValues), 1, 1, 3);
% else
%     try
%         cir=cir.*uint16(~NoValues);
%     catch
%         cir=cir.*uint8(~NoValues);
%     end
% end
% %% Set NoValue cells ==NaN
% 
% % cir(NoValues)=NaN;
% %% compute NDWI or similar index from 3-band image
% [cir_index, initialMask]=waterindex(cir, 'NDWI'); 
% % cir_index=im2uint8(cir_index);
% try
%     cir_index=im2uint8(cir_index+abs(min(min(cir_index))));
% catch
%     cir_index=single(cir_index+abs(min(min(cir_index))));
% end
% 
% %% names
% 
% bw=initialMask;
% mask=NoValues;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
origPixels=sum(sum(bw));
stats=regionprops(mask, 'Area', 'PixelIdxList');
initialstats=regionprops(~bw);
[L, N]=bwlabel(mask);
SE=strel('disk', 1);
% SE=ones(3,3);
for i=1:N
%     disp(i)
%     disp(mean(mean(bw(ring))))
    dil=imdilate(L==i, SE);
    ring=dil & L==0;
    rang=range(bw(ring));
    if rang<=1
%         disp('if')
        stats(i).Fill=true;
        stats(i).FillValue=logical(round(mean(mean(bw(ring)))));
        n=length(stats(i).FillValue);
        bw(stats(i).PixelIdxList)=ones(n, 1)*stats(i).FillValue;
    else
%         disp('else')
        stats(i).FillValue=0;
        
        % unnecessary bc Novalue cells are class as 0
        n=length(stats(i).FillValue);
        bw(stats(i).PixelIdxList)=ones(n, 1)*stats(i).FillValue;
    end
end
filled=bw;
finalPix=sum(sum(filled));
finalstats=regionprops(~filled);
fprintf('\tFilled all NoData cells completely surrounded by water class.\n')
fprintf('\tNumber of cells filled:\t%d\n', finalPix-origPixels)
fprintf('\tNumber of regions filled:\t%i\n', length(initialstats)-length(finalstats))
%% fill



% filled=logical(round(regionfill(uint8(bw), ~mask)));