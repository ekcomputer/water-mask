function complete_region=growUntil(g, spIncl, outputImage, sp_mean, lims)
% (NOT) recursively dilates a superpixel 'image' (graph) for the region containing
%labled sps 'SP_incl' (vector), using graph g (of initial water SPs) until
% condition is met
% lim are limits of growing (fraction) ex: .9, 1.1
% add variance, entropy, or texture image 
% calls SP_dil()


%%
% init
dil_sps=[]; %vector
% completeRegion_old=0;
newring=1.0;
complete_region=spIncl;
% loop
while newring~=0;
    completeRegion_old=complete_region;
    ring=setdiff(SP_dil(g, complete_region), complete_region);
    mean_region=mean(sp_mean(complete_region));
    bounds.a=lims(1)*mean_region;bounds.b=lims(2)*mean_region;
%     fprintf('Region mean= %f\n', mean_region)
%     fprintf('Length of complete_region= %d\n', length(complete_region))
%     fprintf('Length of newring= %d\n', length(newring))
%     fprintf('bounds= %5.0f\t%5.0f\n', bounds(1).a, bounds(1).b)
%      std_region=std(sp_mean(SP_incl)); % single/double prob...

    newring=ring(sp_mean(ring)>=bounds.a & sp_mean(ring)<=bounds.b);
    complete_region=unique([complete_region; newring]);
%     SP_plot(complete_region, stats, size(L_all)); pause(.01)
end