function [regiond,Lnew,L_all] =regionFill(L_all,bw,outputImage, sp_mean, sp_text, cir_index)
% V2.1.5 has different regiond init, and fixes left bar prob...
% V2.1.4 has no shrink in function
% V213 includes mergeRegionsSimple for faster ops
% V2_2_2 includes size limit for growing a small region
% V2_1_1 debugs V2
% V2_1 uses dynamic programming for the averaging
% shrinking is done using global thresh
% region shrinking/growing on superpixels
% L_all is label matrix of SPs, bw is optimal connectivity mask
% outputImage=ndwi or other index (uint8)
% returns regiond, a list of SP corr to classified water, and
% Lnew, a label matrix of conglomerates of SP corr to water
% Elim is Entropy index limit, around 200.

% function calls: edgs2adjList, SP_dil, (SP_plot)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%testing
% L_all=L;
% Elim=f.Elim;
% bounds=f.bounds;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global f
%% shrink based on entropy thresh

%% merge regions to reduce vector sizes
% note: L_all labeling will have gaps...
disp('Merging regions (simple)...')
[L_all, sp_mean, sp_rcount, sp_std, outputImage]=mergeRegions_simple(L_all, bw,cir_index, 'mean');
fprintf('Working with %d superpixels.\n', length(sp_mean))
% sp_rcount=ones(size(sp_mean));
%% continue
waterSP_idx=unique(L_all(bw)); % list of water sp
% waterSP_im=L_all.*double(bw); % label matrix of water sp
L_SP=bwlabel(bw);   % label matrix for all regions in bw
% stats = regionprops(L_in, 'PixelList', 'Centroid',...
%     'Area', 'BoundingBox', 'Perimeter');


%% Convert to graph, extract edge pairs, remove all pairs
% connected to 0 region, shift all values down one to corr to L_all

% g = adjacentRegionsGraph(waterSP_im);
disp('Converting to graph...')
g = adjacentRegionsGraph(L_all);

% indx_toRemove=find(g.Edges.EndNodes==1);
% edges=g.Edges.EndNodes;edges(indx_toRemove,:)=[];edges=edges-1;




%% Loop
disp('Starting loop...')
% imagesc(outputImage)
complete_region=0;
regiond=[]; % grown water body
j=1;
for i=1:max(max(L_SP)) % must ignore single SP regions ?? %change
    fprintf('\nregion %d', i)
    spIncl=unique(L_all(L_SP==i)); % superpixels of given lake region
%     meanIndex=mean(outputImage(L_SP==i));
%     stdIndex=std(double((outputImage(L_SP==i)))); % problem with uint...
%     disp('shrinking...')
%     spIncl=shrinkUntil2(spIncl, sp_mean, sp_entropy, 1.0);
%     disp('Dilating...')
    if 1==1 %length(spIncl)>3*f.sz/f.pArea
        complete_region=growUntil(g, spIncl, outputImage, sp_mean,...
            sp_text, sp_std, sp_rcount, f.bounds);
    else complete_region=spIncl;
    end
%     complete_region=spIncl;
%     ring=setdiff(SP_dil(g, SP_incl), SP_incl);
%     SP_plot(complete_region, stats, size(L_all), 0); pause(.01)
    regiond=unique([complete_region; regiond]);
%     h=subgraph(g, complete_region); 
%     plot(h); 
%     ucc = centrality(h,'betweenness'); histogram(ucc); pause(0.1)
    % for plotting:
% subplot(4,4, mod(j,16)+1); 
% histogram(sp_mean(complete_region), 'BinMethod', 'integers');
% axis([0 256 0 75]); title(num2str(j)); pause(0.1); 
j=j+1;
end
fprintf('\nDone.\n')
figure;
Lnew=SP_plot_raster(regiond, L_all); close
axis image