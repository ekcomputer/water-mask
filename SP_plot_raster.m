function Lnew=SP_plot_raster(SP, L_all, varargin) 
% Lnew=SP_plot_raster(SP, L_all, {comparison, thresh}, 'complete') 
% optional arguments: [string] comparison is either 'lessthan' or 
% 'greaterthan'; [double] thresh is max threshold for viewing non-index SP vectors
% works for SP that includes zeros (mask vector) (filters them out)
% imagesc for lists of superpixel SP (values are indices to L_all,
% given SP label matrix L_all
% returns Lnew, a label matrix of conglomerates of SP corr to water
% only shows binary plots, colorized by SP index
% warning: may change data type/class...
% [idx]=find(SP~=0);
% SP = ismember(SP(idx),unique(L_all)); 
% SP=setdiff(SP,SP(fin));

% inputs (starts w two recurssive branch
if ~isempty(varargin) & strcmp(varargin{1},'lessthan')
    total =[1:max(L_all(:))];
    Lnew=SP_plot_raster(total(SP<varargin{2}), L_all);
elseif ~isempty(varargin) & strcmp(varargin{1},'greaterthan')
    total =[1:max(L_all(:))];
    Lnew=SP_plot_raster(total(SP>varargin{2}), L_all);
else % branch for normal plotting, or complete SP plotting
    % allows for two forms of SP: indices, or mask:
    if range(SP)==1 % if SP is actually a mask
        SP=find(SP);
    end
    pixelIndexList = label2idx(L_all);
    CC=bwconncomp(L_all); % what if some regions are actually adjacent?
    % CC.PixelIdxList{CC.NumObjects}= label2idx(int16(~bw))% puts land regions at end
    if ~isempty(varargin)& strcmp(varargin{1},'complete') % SP is complete list of SP's  
    CC.PixelIdxList=pixelIndexList; % remove zero values ??
    CC.NumObjects=length(CC.PixelIdxList); %+1  
        Lnew=labelmatrix(CC);
        for k = 1:length(CC.PixelIdxList)%max(L_all(:)) < why are  these not same?
            kth_object_idx_list = CC.PixelIdxList{k}; %sum(Lnew(kth_object_idx_list)>0)
            Lnew(kth_object_idx_list) = repmat(double(SP(k)), length(kth_object_idx_list), 1);
%             disp(SP(k))
%             Lnew(kth_object_idx_list) = repmat(2, length(kth_object_idx_list), 1);
        end
    else
        CC.PixelIdxList=pixelIndexList(SP(SP>0)); % remove zero values ??
        CC.NumObjects=length(CC.PixelIdxList); %+1
        Lnew=labelmatrix(CC);
    end
    %% for viewing values
    % for i=1:500 %max(Lnew(:))
    %     n=length(L(L==i));
    %     L(L==i)=repmat(sp_text(i),n,1);
    % end
    %% plot
    axis image
    imagesc(Lnew);
end