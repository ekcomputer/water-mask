% updated to do all computations in the regional script
clearvars -except hl_global hl hl_fused glwd glwd_fused dcs
close all
global env
numperdecade = 5;
regionsStatus=3;
%% Load in the data

load(env.lake_databases,'hl_fused');
load(env.labels_in); 
addpath('power_law_scripts/');

try
    if regionsStatus==2
        for i=1:length(hl_fused)
            hl_fused(i).Region=hl_fused(i).Region2;
        end
    elseif regionsStatus==3
        for i=1:length(hl_fused)
            hl_fused(i).Region=hl_fused(i).Region3_1;
        end
    end
        hl_fused=rmfield(hl_fused, {'Region_1', 'Region2_1', 'Region3_1'});
catch 
    warning('Error.')
end
%% Examine the Fused data for the entire region
disp('Extracting the Fused Dataset');
Fused_area = extractfield(hl_fused,'Area');
% Fused_perim = extractfield(hl_fused,'Perimeter');
Fused_regionnum = extractfield(hl_fused,env.region);
Fused_categorynum = extractfield(hl_fused,env.category);

Regions = sort(unique(Fused_regionnum));
Regions0=[1, env.regions_Q,env.cat_Q]; % Q one
nRegions = length(Regions); 

% [alpha_regional,xmin_regional,pval_regional] = deal(zeros(nRegions,2)); 
ebar_regional = zeros(nRegions,3); 
Fused_regional = cell(nRegions,1); 

for i = Regions0 %1:nRegions
%     if i==21 %yukon flats
%         x1=[yf.Area];
%         p1=[yf.Perimeter];
    if i==1
        Fused_regional{i} = Fused_area;
    elseif i>=22 % shield w and e
        Fused_regional{i} = Fused_area(Fused_categorynum==(i-(min(env.cat_Q)-1)));
    else
        Fused_regional{i} = Fused_area(Fused_regionnum == i); 
    end
%
    
end
labels
%%
Q=Regions0;
for i = Q
    
    [alpha_regional(i,:),xmin_regional(i,:),~,pval_regional(i,:),~,ebar_regional(i,:)] = get_PL_info(Fused_regional{i});
   
end


%% Plotting
close all
set(0,'DefaultLineMarkerSize',6);
set(groot,'defaultAxesFontSize',11); %EK I changed this to 14
set(groot,'defaultLineLineWidth',1);

c=1;
for i = env.regions_Q
%     subplot(2,2,c)
    subplot(4,4,c)
    make_PL_plot(Fused_regional{i},alpha_regional(i,:),xmin_regional(i,:),numperdecade,pval_regional(i,:),ebar_regional(i,:))
    xlabel('Area ($km^2$)'); ylabel('Count'); 
    title(sprintf('%s:\nn = %d, $\\alpha$ = %0.2f',labels{(i)},round(ebar_regional(i,3)), alpha_regional(i,2)));
    set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100 1000])
    c=c+1;
end


pos = [16 8];

%% save
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
saveas(gcf,'Lake_Distributions_Regional.pdf');
saveas(gcf,'Lake_Distributions_Regional.fig');

% save('fitting_data_regional.mat','*_regional');
save('fitting_data_regional.mat','*_regional');
save(env.fit_data_reg,'*_regional');