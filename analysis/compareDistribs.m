% compares GLWD to my results.  input from arcgis

% TODO: figure out how not to distort plots- combine datasets before
% calling ecdf- magnitude matters whne rescaling log-log plots

%% run previous script to get DCS data
analyzeWaterDistribution

%% set params

noDuplicates=0; % don't plot same lakes on AS and GLWD curves
noCurveOverlap=1; % cut off high end of airswot
useLog=0; % log transform data before plotting
region=1; % use 1 for all
%% load GLWD data
g1_in='F:\GLWD\out\GLWD_1_lent_intrsct.shp';
g2_in='F:\GLWD\out\GLWD_2_lent_intrsct.shp';

g1x_in='F:\GLWD\out\GLWD_1_lent_xbounds.shp';
g2x_in='F:\GLWD\out\GLWD_2_lent_xbounds.shp';

[~, g1a]=shaperead(g1_in);
[~, g2a]=shaperead(g2_in);

[~, g1xa]=shaperead(g1x_in);
[~, g2xa]=shaperead(g2x_in);

%% preprocess glwd data
f.ar1=[g1a.AREA_SKM];
f.per1=[g1a.PERIM_KM];
f.ar2=[g2a.AREA_SKM];
f.per2=[g2a.PERIM_KM];

f.arx1=[g1xa.AREA_SKM];
f.perx1=[g1xa.PERIM_KM];
f.arx2=[g2xa.AREA_SKM];
f.perx2=[g2xa.PERIM_KM];

glwd.ar=[f.ar1, f.ar2];
glwd.per=[f.per1, f.per2];

glwd.xar=[f.arx1, f.arx2];
glwd.xper=[f.perx1, f.perx2];
clear f
%% make plots
    % cumulative area
    figure; hold on
% [N,edg]=histcounts(abun_rshp.ar(rshp_msk)/1e6, ev_ar, ...
% 'Normalization', 'cumcount');
% plot(edg(1:end-1), max(N)-N); xlabel('Area ($km^2$)'); ylabel('Count of lakes greater than given area');

    % prep plotting data
    
if ~noDuplicates
    glwd.use=glwd.ar;
else
    glwd.use=glwd.xar; % use only cross boundaries area
end

dcs_area=abun_rshp.ar(rshp_msk(:,region))/1e6;
if noCurveOverlap
    dcs_area=dcs_area(dcs_area<min(glwd.use));
end

if useLog;
   glwd.use=log10(glwd.use+1e-14);
   dcs_area=log10(dcs_area+1e-14);
end

    % plot pareto fit

dcs_area=[dcs_area, glwd.ar];
h.bnd=min(glwd.ar); % boundary between two plots
hold off
[CdfY,CdfX] = ecdf(dcs_area,'Function','survivor'); 
plot(CdfX(CdfX<h.bnd), CdfY(CdfX<h.bnd), '--b', 'LineWidth', 2.5)
hold on
plot(CdfX(CdfX>=h.bnd), CdfY(CdfX>=h.bnd), ':k') %rescale(CdfY_G, min(CdfY_G), 0.01961))
legend({'AirSWOT camera lakes', 'GLWD lakes'}, 'location', 'best')

hold off
if useLog
    xlabel('$Log_{10}$(Area) ($km^2$)'); ylabel('$Log_{10}$(Number of lakes of greater area)');
else
    xlabel('Area ($km^2$)'); ylabel('Number of lakes of greater area');
end
title('Lake area distributions at two scales')
if ~useLog
    
    set(gca, 'YScale', 'log', 'XScale', 'log')
else
    set(gca, 'YTickLabel', [])
end