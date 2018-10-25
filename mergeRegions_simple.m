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

l.before=max(L_all(:));
L_all=~bw.*L_all+(max(L_all(:))+1)*bw+bwlabel(bw); % leaves gaps
L_all=fillLabelGaps(L_all); % fills gaps
% L_all=relabel(L_all); % relabels in viert-horz order (matlab style)
% stats=regionprops(L_all, 'area', 'pixelidxlist');
stats=regionprops(L_all, 'area');
keeperIndices=find([stats.Area]~=0);
stats=stats(keeperIndices);
l.len=length(stats); clear stats
%% loop

outputImage_out = zeros(size(cir_index),'like',cir_index);
sp_out=zeros(l.len,1,'like',cir_index); %sp_dev=zeros(N,1,'like',A);
idx = label2idx(L_all);
for i = 1:length(idx) % < here!
    cir_index_Idx = idx{i};
    sp_rcount(i)=length(cir_index(cir_index_Idx)); %count of superpixel (in pixels)
    sp_std(i)=std(double(cir_index(cir_index_Idx)), 'OmitNaN'); %std dev of superpixel (in pixels)
    if isempty(varargin) | strcmp(varargin{1}, 'mean')
        outputImage_out(cir_index_Idx) = mean(cir_index(cir_index_Idx), 'OmitNaN');
        sp_out(i)=mean(cir_index(cir_index_Idx), 'OmitNaN'); %mean value of superpixel
%         mean(cir_index(cir_index_Idx));
%         i;
%         imagesc(outputImage_out);
    elseif strcmp(varargin{1}, 'std')
        outputImage_out(cir_index_Idx) = std(double(cir_index(cir_index_Idx)), 'OmitNaN');
        sp_out(i)=std(double(cir_index(cir_index_Idx)), 'OmitNaN'); 
    elseif strcmp(varargin{1}, 'count')
        outputImage_out(cir_index_Idx) = sum(sum(double(cir_index(cir_index_Idx))), 'OmitNaN');
        sp_out(i)=sum(sum(double(cir_index(cir_index_Idx)), 'OmitNaN')); 
    end
end   
sp_rcount=sp_rcount';
l.after=length(unique(L_all));
disp('Done merge Regions Simple.')
fprintf('\tMerged %d regions.\n',l.before-l.after)