% to run after waterDistribution.m
% fixes entries with no water - just use once!

clear
struct_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\distrib.mat';
load(struct_out)

%%
for i=1:length(abun)
    fprintf('%d   %d\n', i,length(abun(i).stats) )
    if isempty(abun(i).stats)
        abun(i).stats=rmfield(abun(i).stats, 'PixelList');
        abun(i).stats(1).Area=[];
        abun(i).stats(1).Centroid=[];
        abun(i).stats(1).MajorAxisLength=[];
        abun(i).stats(1).Perimeter=[];
        abun(i).stats(1).lat=[];
        abun(i).stats(1).long=[];
        abun(i).stats(1).minBoundRadius=[];
        abun(i).stats(1).LehnerDevel=[];
        abun(i).stats(1).SDF=[];
    end
end

%% save
save(struct_out, 'abun')
fprintf('Saving to %s\n', struct_out)
disp('Done.')