function dil_sps=SP_dil(g, SP_incl)
% dil_sps=SP_dil(g, SP_incl)
% "Superpixel dilation": dilates a superpixel 'image' (graph) for the region containing labled
% sps 'SP_incl' (vector), using graph g (of initial water SPs)

dil_sps=[]; %vector
for i=1:length(SP_incl)
%     dil_sps=unique([dil_sps; neighbors(g,SP_incl(i))]);
    dil_sps=unique([SP_incl; dil_sps; neighbors(g,SP_incl(i))]);

end

%% Plot check
% 
% for i= 1:length(dil_sps)
%     axis([0 500 750 1250])
%     imagesc(L_all==dil_sps(i))
% end