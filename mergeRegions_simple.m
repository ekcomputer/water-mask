function [L_all, sp_out, sp_rcount, sp_std, outputImage_out]=mergeRegions_simple(L_all, bw,...
    cir_index, varargin)
% V3.3.2 fixes memory usage bug introduced by fixing left bar problem!
% V3.3 fixes left bar problem
% V32 fills gaps in label matrix, L, labeling
% V3_1 is debugged V3, but still leaves gaps in labeling (note:
% outputImage for 'count' doesn't work
% V3 has simplified and updated output to include SP_regioncount_out
% merges regions in L_all (label matrix) based on regions in bw (binary
% image).  Returns simplified SP image (outputImage), vector SP_out, 
% SP_rcount (which is vector of number of combined SPs), and
% updated label matrix (L_all_out)
% varargin can be 'mean' (default), 'std' or 'count'
% calls fillLabelGaps

    % first, create unique index for new L values
l.before=max(L_all(:));
% L_all=uint32(~bw).*L_all+(max(L_all(:))+1)*uint32(bw+bwlabel(bw)); % leaves gaps
L_all=uint32(~bw).*L_all+(max(L_all(:))+1)*ones(size(L_all), 'like', L_all).*uint32(bw)+...
    uint32(bw+bwlabel(bw)); % leaves less gaps
L_all=fillLabelGaps(L_all); % fills gaps
% L_all=relabel(L_all); % relabels in viert-horz order (matlab style)
% stats=regionprops(L_all, 'area', 'pixelidxlist');
stats=regionprops(L_all, 'area');
keeperIndices=find([stats.Area]~=0);
stats=stats(keeperIndices);
l.len=length(stats); clear stats
%% loop
    % init vars
% outputImage_out = zeros(size(cir_index),'like',cir_index);
% sp_out=zeros(l.len,1,'like',cir_index); %sp_dev=zeros(N,1,'like',A);
% sp_rcount=zeros(l.len,1);
% sp_std=zeros(l.len,1);
    % loop

idx = label2idx(L_all);
cellmean=@(x)mean(cir_index(x), 'OmitNaN');
sp_out=uint8(cellfun(cellmean, idx, 'UniformOutput', true));
cellstd=@(x)std(double(cir_index(x)), 'omitnan');
sp_std=cellfun(cellstd, idx, 'UniformOutput', true);
cellcount=@(x)numel(~isnan(cir_index(x)));
sp_rcount=cellfun(cellcount, idx, 'UniformOutput', true);
outputImage_out=SP_plot_raster(sp_out, L_all, 'complete', 'noplot');

sp_rcount=sp_rcount';
sp_out=sp_out';
sp_std=sp_std';
l.after=length(unique(L_all));
disp('Done merge Regions Simple.')
fprintf('\tMerged %d regions.\n',l.before-l.after)