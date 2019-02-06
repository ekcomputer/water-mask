% compares GLWD to my results.  input from arcgis

% TODO: figure out how not to distort plots- combine datasets before
% calling ecdf- magnitude matters whne rescaling log-log plots

%% run previous script to get DCS data
analyzeWaterDistribution
fprintf('Note: Max DCS lake size: %d\n', maxSize)
%% set params

noDuplicates=1; % don't plot same lakes on AS and GLWD curves
noCurveOverlap=1; % cut off high end of airswot
useLog=0; % log transform data before plotting
region=1; % use 1 for all
fusedPlot=1; % make seamless plot with redundancy
overlapPlots=0; % plot both on same axis % dont use if fusefPlot=1
overlapPlots2=0; % alternate way
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
dcs_perim=abun_rshp.per(rshp_msk(:,region))/1e3;
% if noCurveOverlap
%     dcs_area=dcs_area(dcs_area<min(glwd.use));
% end

if useLog;
   glwd.use=log10(glwd.use+1e-14);
   dcs_area=log10(dcs_area+1e-14);
end




    % plot pareto fit
fused_area=[dcs_area, glwd.use];
    % create mask- 1 means dcs obs, 0 means glwd obs
h.fuseMask=false(size(fused_area)); h.fuseMask(1:length(dcs_area))=true;
if fusedPlot % make seamless plot with possible redundancy
    
    if noCurveOverlap
                % premask with bound
        h.bnd=min(glwd.ar); % boundary between two plots
        h.bnd=0.5; % not used
        h.low=0.1; % low bound
        h.high=1;
        dcs_area_dcs=dcs_area(dcs_area<h.bnd);
        glwd.use=glwd.use(glwd.use>=h.bnd);
    else dcs_area_dcs=dcs_area;
    end
    hold off
    [CdfY,CdfX,h.l,h.h] = ecdf(fused_area,'Function','survivor'); 
    plot(CdfX(CdfX<h.low), CdfY(CdfX<h.low), '--b', 'LineWidth', 2.5)
    hold on
    plot(CdfX(CdfX>=h.high), CdfY(CdfX>=h.high), ':k') %rescale(CdfY_G, min(CdfY_G), 0.01961))
    plot(CdfX(CdfX<h.high & CdfX>=h.low), CdfY(CdfX<h.high & CdfX>=h.low), '-.', 'color',[0.5 0.5 0.5]')
    plot(CdfX, h.l, ':k', 'LineWidth', 1); plot(CdfX, h.h, ':k','LineWidth', 1);
    legend({'AirSWOT camera lakes', 'GLWD lakes', 'mixed', '95 \% confidence'}, 'location', 'best')
    hold off
elseif overlapPlots
%     [CdfY,CdfX] = ecdf(dcs_area,'Function','survivor'); 
%     [CdfY_G,CdfX_G] = ecdf(glwd.use,'Function','survivor');
    [CdfY,CdfX] = ecdf(fused_area,'Function','survivor', 'censoring', ~h.fuseMask);
    [CdfY_G,CdfX_G] = ecdf(fused_area,'Function','survivor', 'censoring', h.fuseMask);
    plot(CdfX, CdfY, '--b', 'LineWidth', 2.5)
    hold on
    plot(CdfX_G(CdfY_G<1), CdfY_G(CdfY_G<1), ':k') %rescale(CdfY_G, min(CdfY_G), 0.01961))  
    hold off
    
elseif overlapPlots2
    ev_ar=[0:0.0001:max(dcs_area)]; % redefine using linear bins
    ev_ar_G=[min(glwd.use):0.1:max(glwd.use)];
    [N,edg]=histcounts(dcs_area, ev_ar, 'Normalization', 'cumcount');
    [N_G,edg_G]=histcounts(glwd.use, ev_ar_G, 'Normalization', 'cumcount');
%     for n=1:length(ev_ar)
%        edg(n)=ev_ar(n);
%        N(n)=sum(dcs_area>edg(n));
%     end
    plot(edg(1:end-1), max(N)-N, '--b'); xlabel('Area ($km^2$)'); ylabel('Count of lakes greater than given area');
    hold on
    plot(edg_G(1:end-1), max(N_G)-N_G, ':r');
    legend({'AirSWOT camera lakes', 'GLWD lakes'}, 'location', 'best')
    hold off
end

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

%% re-do calcs using fused GLWD dataset
    i=1;
    fusedTotal.pers=[dcs_perim, glwd.xper];
    fusedTotal.ars=[dcs_area, glwd.xar];
    fusedTotal.lim=total.lim;
    fusedTotal.totalLenticAr=sum(fusedTotal.ars);
    fusedTotal.totalLenticPer=sum(fusedTotal.pers);
    fusedTotal.minSize=total(i).minSize;
    fusedTotal.maxSize=total(i).maxSize;  
    fusedTotal.count=length(fusedTotal.ars);
    fusedTotal.perUnder01=sum(fusedTotal.ars<0.01)/fusedTotal.count;
    fusedTotal.perUnder001=sum(fusedTotal.ars<0.001)/fusedTotal.count;
    fusedTotal.perUnder0001=sum(fusedTotal.ars<0.0001)/fusedTotal.count;
    
    fusedTotal.ArPerUnder01=sum(fusedTotal.ars(fusedTotal.ars<0.01))/fusedTotal.totalLenticAr;
    fusedTotal.ArPerUnder001=sum(fusedTotal.ars(fusedTotal.ars<0.001))/fusedTotal.totalLenticAr;
    fusedTotal.arPerUnder0001=sum(fusedTotal.ars(fusedTotal.ars<0.0001))/fusedTotal.totalLenticAr;
    
    fusedTotal.PerPerUnder01=sum(fusedTotal.pers(fusedTotal.ars<0.01))/fusedTotal.totalLenticPer;
    fusedTotal.PerPerUnder001=sum(fusedTotal.pers(fusedTotal.ars<0.001))/fusedTotal.totalLenticPer;
    fusedTotal.PerPerUnder0001=sum(fusedTotal.pers(fusedTotal.ars<0.0001))/fusedTotal.totalLenticPer;
  