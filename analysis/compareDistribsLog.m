% compares GLWD to my results.  input from arcgis

%% run previous script to get DCS data
analyzeWaterDistribution

%% load GLWD data
g1_in='F:\GLWD\out\GLWD_1_lent_intrsct.shp';
g2_in='F:\GLWD\out\GLWD_2_lent_intrsct.shp';

[~, g1a]=shaperead(g1_in);
[~, g2a]=shaperead(g2_in);

%% preprocess glwd data
f.ar1=[g1a.AREA_SKM];
f.per1=[g1a.PERIM_KM];
f.ar2=[g2a.AREA_SKM];
f.per2=[g2a.PERIM_KM];

glwd.ar=[f.ar1, f.ar2];
glwd.per=[f.per1, f.per2];
clear f
%% make plots
    % cumulative area
    figure; hold on
[N,edg]=histcounts(abun_rshp.ar(rshp_msk)/1e6, ev_ar, ...
'Normalization', 'cumcount');
plot(edg(1:end-1), max(N)-N); xlabel('Area ($km^2$)'); ylabel('Count of lakes greater than given area');


%         figure(2); hold on



% plot pareto fit
if 1==1
    hold off
    [CdfY,CdfX] = ecdf(log10(abun_rshp.ar(rshp_msk)/1e6+1e-8),'Function','survivor'); 
    plot(CdfX, CdfY)
    hold on
    [CdfY_G,CdfX_G] = ecdf(log10(glwd.ar+1e-8),'Function','survivor'); 
    plot(CdfX_G, CdfY_G) %rescale(CdfY_G, min(CdfY_G), 0.01961))
        % fits
    pfit{1}=polyfit(CdfX, CdfY, 1)
    pfit{2}=polyfit(CdfX_G, CdfY_G, 1)
    line_x{1}=@(x) x*pfit{1}(1) + pfit{1}(2);
    line_x{2}=@(x) x*pfit{2}(1) + pfit{2}(2);
    fplot(line_x{1})
    fplot(line_x{2})
    legend({'AirSWOT camera lakes', 'GLWD lakes'}, 'location', 'best')
end
hold off
xlabel('Area ($km^2$)'); ylabel('Number of lakes of greater area');
title('Lake area distributions at two scales')
% set(gca, 'YScale', 'log', 'XScale', 'log')
if i==i_end % last time
    if plotRegions
        legend(labels, 'location', 'best', 'FontSize', 15)
    else
    end
end