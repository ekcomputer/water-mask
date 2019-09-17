% function for loading and plotting area change stats

clear
matPath_base='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\areaChange\mat\';
matPath=[matPath_base, 'Change.mat'];
load(matPath); % gives me stats_all structure


%%% choose params
id=57; % file id, if no loop
minSize=150; % min water body size to consider in plot
%%%

%%
pchange_all=[];
change_all=[];
mask_all=[];
mask_all=[];
for n=1:length(stats_all)
% for n=id
    if isempty(stats_all(n).change)
        continue
    end
    pchange=[stats_all(n).change.pchange];
    change=[stats_all(n).change.change];
    size_a=[stats_all(n).change.size_a];
    size_b=[stats_all(n).change.size_b];

    mask=pchange>-9999 & pchange< 9999 & size_a>minSize & size_b>minSize;
    histogram(pchange(mask), 'BinMethod', 'auto')
    title(stats_all(n).change_sum.file_b(18:26))
%     pause(0.5) % drawnow
    pchange_all=[pchange_all,pchange]; %concat all
    change_all=[change_all,change]; %concat all
    mask_all=[mask_all, mask]; %concat all
end
figure
histogram(pchange_all(logical(mask_all)), 'BinWidth', 20, 'BinLimits', [-100 1000] )
title('Percent Change (all)')
xlabel('\%')
ylabel('Count')

    % absolute change
figure
histogram(change_all(logical(mask_all)), 'BinWidth', 20, 'BinLimits', [-100 1000] )
title('Abolute Change (all)')
xlabel('$Meters^2$')
ylabel('Count')

%% reshape struct
change_all=vertcat(stats_all.change);
mask=[change_all.pchange]>-9999 & [change_all.pchange]< 9999 &...
    [change_all.size_a]>minSize & [change_all.size_b]>minSize;
figure
histogram([change_all(mask).size_a]); set(gca,'YScale','log')
xlabel('Area')
ylabel('Count')

figure
subplot(211)
plot([change_all(mask).lat], [change_all(mask).size_a], '.');
xlabel('Latitude')
ylabel('Area')

subplot(212)
plot([change_all(mask).long], [change_all(mask).size_a], '.');
xlabel('Longitude')
ylabel('Area')