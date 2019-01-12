% script to extract values from vic or reannalysis data and compare to WSA
% change within that grid cell
clear
close all
% f.shp_DCS_loc='J:\output\dcs_tileIndex\DCS_Index_simpl.shp';
% [shp, att]=shaperead(f.shp_DCS_loc);

struct_in='J:\output\analysis\distrib.mat';
pe_in='J:\AboveP-E\2017_geotiff_WGS84_EPSG4326\ABoVE.P-E.7_2017.tif';
pairs_path='J:\Final\analysis\pairs.mat';
change_pth='J:\Final\analysis\areaChange\shp\aggregate\water_change_2.shp';
load(struct_in); load(pairs_path);
[pe, R]=geotiffread(pe_in);
[shp, att]=shaperead(change_pth);

minSize=4000;
imagesc(flipud(pe))
pe=flipud(pe); % necessary?
%% loop
for i=1:length(att)
    [x,y]=geographicToIntrinsic(R, att(i).lat, att(i).long);
    mm(i)=pe(round(y),round(x));
end

%% filter out appears and disappears (pre-process)
pchange=[att.pchange3];
pchange_no0=pchange(pchange<100 & pchange>-100); mm_no0=mm((pchange<100 & pchange>-100));
pchange_mskd=pchange([att.size_a]>minSize | [att.size_b]>minSize); mm_mskd=mm([att.size_a]>minSize | [att.size_b]>minSize);

%% plot

figure(1)
plot(mm, pchange, '.')
xlabel('P-E (mm)'); ylabel('Percent change (\%)')
title('No filtering')

figure(2)
plot(mm_no0, pchange_no0, '.')
xlabel('P-E (mm)'); ylabel('Percent change (\%)')
title('Filtering: only lakes that persist')

figure(3)
plot(mm_mskd, pchange_mskd, '.')
xlabel('P-E (mm)'); ylabel('Percent change (\%)')
title(['Only lakes above ', num2str(minSize/1000000), '$km^2$'])