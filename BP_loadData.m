function [cir_index, NoValues,waterFlag, medWaterIndex]= BP_loadData(cir, varargin)
% uses median instead of mean for water index... using quanitles 
% varargin can be index type as first, then saturation percent for index as
% second
% rescales histogram for cir_index ever-so-slightly
% options:  case 'NDWI''BNDWI''MNDWI''NDVI''WI''MAXNDWI''SPAN''MAXWI'
% 'MINNDWI''MINWI''BNDWI_3_100''BNDWI_1_10'  'AWEI' 'IR'
% calls waterindex
% i is training file number
% varargin = {water_index}
% addpath D:\Dropbox\Matlab\Above\
waterFlag=[-99 -99 -99];
if length(varargin)>0
    waterIndex=varargin{1}
else waterIndex='NDWI'
end
if length(varargin) >1 && strcmp(varargin{2}, 'satPercent')
    Tol=varargin{3};
else
    Tol=0.00005;
end
% Rewritten to detect SP on masked image
% Next: gradient to detect shoulders
% function labeled=super_1.m(cir, waterIndex)
% script to use OBIA and superpixels to classify open water extent

%%%%%%%%%%%%%%%%%%%%non- function params
% clear
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set NoValue cells ==0
% NoValues=or(cir(:,:,1)==192 & cir(:,:,2)==192 & cir(:,:,3)==192,...
%     cir(:,:,1)==1 & cir(:,:,2)==1 & cir(:,:,3)==1); %255
data_type=class(cir);
if strcmp(data_type, 'int16')
    NoValues=cir(:,:,1)==-9999 & cir(:,:,2)==-9999 & cir(:,:,3)==-9999;
    cir=im2uint16(rescale(cir, 'InputMin',0, 'InputMax',32767));
else
    NoValues=or(cir(:,:,1)==65535 & cir(:,:,2)==65535 & cir(:,:,3)==65535,...
        cir(:,:,1)==0 & cir(:,:,2)==0 & cir(:,:,3)==0); %255
end
if size(cir, 1)*size(cir, 2)==sum(NoValues(:)) % condition for novalue tile
    cir_index=zeros(size(NoValues));
    waterFlag(1)=0; medWaterIndex=-99;
    disp('NoData tile.')
    return
end
b1=cir(:,:,1);
b2=cir(:,:,2);
b3=cir(:,:,3);

b1(NoValues)=0;
b2(NoValues)=0;
b3(NoValues)=0;

cir=cat(3, b1, b2, b3);
clear b1 b2 b3
% if strcmp(class(cir),'single')
%    cir=cir.*repmat(single(~NoValues), 1, 1, 3);
% else
%     try
%         cir=cir.*uint16(~NoValues);
%     catch
%         cir=cir.*uint8(~NoValues);
%     end
% end
%% Set NoValue cells ==NaN

cir(NoValues)=NaN;
%% compute NDWI or similar index from 3-band image
global f
[cir_index]=waterindex(cir, waterIndex, NoValues); 
f.origLevel=graythresh(cir_index(~NoValues))
if f.plot
    figure; imagesc(imoverlay(cir, cir_index>f.origLevel)); axis image
    title('NDWI threshold from raw image')
end
%% Condition for no water
medWaterIndex=median(cir_index(:), 'omitnan')
% waterFlag(2)=median(cir_index(cir_index>quantile(cir_index(:), 0.9999))); % median "water"
waterFlag(2)=median(maxk(cir_index(:), f.minAreaFact*f.minSize/f.pArea)); % median "water"
waterFlag(3)=median(mink(cir_index(:), f.minAreaFact*f.minSize/f.pArea)); % median "land"
if waterFlag(2)< f.NDWIWaterAmount % comparing median water %no water
    waterFlag(1)=0;
    disp('No water detected.')
elseif waterFlag(3)> f.NDWILandAmount % comparing median land % no land
    waterFlag(1)=2;
    disp('No land detected.')
else
    waterFlag(1)=1; % there is water
end
if f.plot
    figure; subplot(311);
    histogram(cir_index(~NoValues))
    hold on; plot([f.NDWILandAmount f.NDWIWaterAmount], [0, 0], 'bV'); hold off
    title('Raw histogram')
end
% subplot(122)
% imagesc(cir)


%% positive-valued and rescaled assume NoData value is same in all bands
cir_index_pos=rescale(cir_index+abs(min(0, min(cir_index(:))))); 

% enhanced to saturate top n % of pixels
cir_index_enh=imadjust(cir_index_pos, stretchlim(cir_index_pos(~NoValues), Tol)  );
if f.plot
    subplot(312)
    histogram(cir_index_enh(~NoValues))
    title('Stretched .0005\%')
end
%% Sigmoid transform on image
    % problem is that transform can remove 'bumps' in histogram, which
    % messes up opotmizeConn
% cir_index_enh=ParabAdjust(cir_index_enh, NoValues, 1.3733, -2.4,  3.7333, -1.7067);

%% convert to uint8
switch data_type % uncessary if then tree...
    case 'uint8'
        cir_index=im2uint8(cir_index_enh);
    case 'uint16'
        cir_index=im2uint8(cir_index_enh);
    case 'int16'
        cir_index=im2uint8(cir_index_enh);
    case 'single'
%         warning('EK error')
        cir_index=im2uint8(cir_index_enh);
%         cir_index=(cir_index+abs(min(min(cir_index)))); %?
    case 'double'
        error('error')
%         cir_index=(cir_index+abs(min(min(cir_index)))); %?
end
if f.plot
    subplot(313)
    histogram(cir_index(~NoValues))
    title('Scaled to int8')
end
% figure
%% binarize and visualize
% level=graythresh(cir_index(~NoValues));
% initialMask=imbinarize(cir_index, level-.005);
% ol=imoverlay(cir, initialMask);
% figure
% imagesc(ol); axis image;title('Initial Mask from generously binarizing index')