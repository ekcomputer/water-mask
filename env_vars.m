% environment variables
% specify most recent versions of files
global env
try
    env.regions_Q=2:14; % for power law plots and table
    env.cat_Q=22:26;%[1:15]; % labels(Q) to plot

    % file paths
    env.saved_dir='D:\GoogleDrive\Research\Lake distributions\savedData\';
    env.tbl_out=[env.saved_dir, 'out\LakeMorphology_4_1.xlsx'];
    env.geom_in=[env.saved_dir, 'Geom2.mat'];
    env.labels_in='D:\GoogleDrive\Research\Lake distributions\regionLabels4.mat';
    env.labels_exp_in='D:\GoogleDrive\Research\Lake distributions\regionLabels4_expl.mat';
    env.shp_in='F:\AboveDCSRasterManagement\DataMaskShapes\Regions\DCS_regions_4.shp';
    env.shp_diss_in='F:\AboveDCSRasterManagement\DataMaskShapes\Regions\DCS_regions_4_Dissolve.shp';
    env.wc_complete='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_combined\WC_complete_water_diss.shp';
    env.fused='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\dcs_fused_hydroLakes.shp';
    env.fused_xs='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\dcs_fused_hydroLakes_xs.shp'; % fused lakes intersected with study region
    env.lake_databases='D:\GoogleDrive\Research\Lake distributions\LakeDatabases.mat';
    env.analyzeWaterDistribution=[env.saved_dir, 'analyzeWaterDistribution_4_1.mat'];
    env.fit_data=[env.saved_dir, 'fitting_data_det2.mat'];
    env.fit_data_reg=[env.saved_dir, 'fitting_data_regional_4_det.mat'];
    env.plotGCPAcc=[env.saved_dir, 'plotGCPAcc.mat'];
    env.region='Region4';
    env.category='Category4';
    env.plinfo.repnum=1000;
    env.plinfo.samp_num=1000;
    env.plinfo.semi_param_number=1000;
    env.buffer_dir='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\summ_output';  % buffer analysis directory
    env.hl_fused_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\summ_output\dcs_fused_hydroLakes_buf_10_sum.shp';
    env.bulk_plot_dir='D:\pic\bulk\';

catch
    warning('Trouble loading env. variables (Ethan)')
end