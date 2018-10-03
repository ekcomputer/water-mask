function complete_region=growUntil(g, spIncl, outputImage, sp_mean, sp_rcount, lims)
% V5.1 fixes a bug in calculating region mean vs region index...
% V5 uses proper weighting.  .count refers to pixels
% V4 makes uses of simply-merged regions, so needs a region size input
% V2 use alternative to unique() function and fixes mean bug
% V1 uses dynamic programming
% (NOT) recursively dilates a superpixel 'image' (graph) for the region containing
%labled sps 'SP_incl' (vector), using graph g (of initial water SPs) until
% condition is met
% lim are limits of growing (fraction) ex: .9, 1.1
% note spIncl gives indexes to sp_mean;
% add variance, entropy, or texture image 
% calls SP_dil() and fastSetdiff()


%%
% init
dil_sps=[]; %vector
% completeRegion_old=0;
newring.idxs=1.0;
complete_region=spIncl; clear spIncl % units of index!
firstTime=true;
% loop
while newring.idxs~=0;
    ring.idxs=setdiff(SP_dil(g, complete_region), complete_region);
%     ring.idxs=otherFastSetdiff(SP_dil(g, complete_region), complete_region);
    region.idxs=complete_region; %
    if firstTime
        region.mean=mean(sp_mean(region.idxs)); % no weighting here
    else
        region.mean=(region.count*region.mean+newring.count*newring.mean)/(region.count+newring.count);
    end
    region.count=length(region.idxs)*sum(sp_rcount(region.idxs));
    bounds.a=lims(1)*region.mean  ;  bounds.b=lims(2)*region.mean;
%     fprintf('Region mean= %f\n', mean_region)
%     fprintf('Length of complete_region= %d\n', length(complete_region))
%     fprintf('Length of newring= %d\n', length(newring))
%     fprintf('bounds= %5.0f\t%5.0f\n', bounds(1).a, bounds(1).b)
%      std_region=std(sp_mean(SP_incl)); % single/double prob...

    newring.idxs=ring.idxs(sp_mean(ring.idxs)>=bounds.a & sp_mean(ring.idxs)<=bounds.b);
    newring.count=length(newring.idxs)*sum(sp_rcount(newring.idxs));
    newring.mean=mean(double(sp_mean(newring.idxs)).*sp_rcount(newring.idxs));
    complete_region=unique([complete_region; newring.idxs]);
    if newring.count>0 % &~firstTime
        fprintf(' %d', length(newring.idxs))
    end
    firstTime=false;
%     SP_plot(complete_region, stats, size(L_all)); pause(.01)
end
% complete_region=unique(complete_region);