clearvars -except hl_global hl hl_fused glwd glwd_fused dcs
close all
numperdecade = 5; 
%% Load in the data

load('D:\GoogleDrive\Research\Lake distributions\LakeDatabases.mat','hl_fused');
addpath('power_law_scripts/');

%% Look at the HydroLakes Dataset

disp('Extracting the HL Dataset'); 
HL_area = extractfield(hl_global,'Lake_area');
[HL_alpha,HL_xmin,~,HL_pval,~,HL_ebar] = get_PL_info(HL_area); 

%% Examine the Hydrolakes Dataset for the Obs region
disp('Extracting the Local Hydrolakes Dataset'); 
HL_local_area = extractfield(hl,'Lake_area');
[HL_local_alpha,HL_local_xmin,~,HL_local_pval,~,HL_local_ebar] = get_PL_info(HL_local_area); 

%% Examine the DCS data for this region
disp('Extracting the DCS Dataset'); 
DCS_area = extractfield(dcs,'Area');
[DCS_alpha,DCS_xmin,~,DCS_pval,~,DCS_ebar] = get_PL_info(DCS_area); 

%% Examine the Fused data for the entire region
disp('Extracting the Fused Dataset'); 
Fused_area = extractfield(hl_fused,'Area');
Fused_perim = extractfield(hl_fused,'Perimeter'); 
[Fused_alpha,Fused_xmin,~,Fused_pval,~,Fused_ebar] = get_PL_info(Fused_area); 


%% Plotting
close all

subplot(231)
make_PL_plot(HL_area,HL_alpha,HL_xmin,numperdecade,HL_pval,HL_ebar)
xlabel('Area ($km^2$)'); ylabel('Count'); title(sprintf('All Hydrolake Data: n = %d',round(HL_ebar(3))));
subplot(232)
make_PL_plot(HL_local_area,HL_local_alpha,HL_local_xmin,numperdecade,HL_local_pval,HL_local_ebar)
xlabel('Area ($km^2$)'); ylabel('Count'); title(sprintf('Local Hydrolake Data: n = %d',round(HL_local_ebar(3))));
subplot(233)
make_PL_plot(DCS_area,DCS_alpha,DCS_xmin,numperdecade,DCS_pval,DCS_ebar)
xlabel('Area ($km^2$)'); ylabel('Count'); title(sprintf('All DCS Data: n = %d',round(DCS_ebar(3))));
subplot(212)
make_PL_plot(Fused_area,Fused_alpha,Fused_xmin,numperdecade,Fused_pval,Fused_ebar)
xlabel('Area ($km^2$)'); ylabel('Count'); title(sprintf('%s:\nn = %d, $\\alpha$ = %0.2f','Fused water bodies',round(Fused_ebar(3)), Fused_alpha(2)));
% title(sprintf('Fused Data: n = %d',round(Fused_ebar(3))));
set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100 1000])

pos = [12 8]; 
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
% saveas(gcf,'Lake_Distributions.pdf'); 
% 
% save('fitting_data.mat','*_xmin','*_alpha','*_area','*_pval','*_ebar'); 
save('D:\GoogleDrive\Research\Lake distributions\savedData\fitting_data.mat','*_xmin','*_alpha','*_area','*_pval','*_ebar');