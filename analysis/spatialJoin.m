% trying to properly do the buffer analysis which arc couldn't do

clear; close all
file_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\iterations\dcs_fused_hydroLakes_buf_50_sum.shp';
file_inF='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\dcs_fused_hydroLakes.shp';

[shp,att]=shaperead(file_in);
[shpF,attF]=shaperead(file_inF); % fused

%% convert 

pshpF=polyshape([shpF.X], [shpF.Y]); % slooooooow
%% loop
for i=1:length(shp)
    pshp(i)=polyshape(shp(i).X, shp(i).Y); % slow
    disp(i)
    ar(i).Area=area(intersect(pshpF,pshp(i))); % area in m
    ar(i).Region4=att(i).Region4;
    ar(i).Category4=att(i).Category4;
end


%% save
save('D:\GoogleDrive\Research\Lake distributions\WaterBodies.mat', 'ar')
disp('File saved.')
%% TODO

% load new file_in     >  run loop overnight
