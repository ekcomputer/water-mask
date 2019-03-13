% script to compute areas of each region and water fractions
clear 
global env
disp('Loading data...')
shp_in=env.shp_in;
load(env.labels_in);
saved_dir=env.saved_dir;
% [S,attS]=shaperead(shp_in);
[SD,attD]=shaperead(env.shp_diss_in); % dissolved shapefile
[~,attWCC]=shaperead(env.wc_complete); % all water (wc_complete_water)
[~,attF]=shaperead(env.fused_xs); % all lakes (with extra areas clipped)
Q=[1,env.regions_Q, env.cat_Q];
disp('Data loaded.')
%% loop for areas
for i= 1:length(SD)%Q
   SD(i).shp=polyshape(SD(i).X, SD(i).Y) 
   SD(i).area=area(SD(i).shp)/1e6;
end

for i=Q
    if i==1
        geom(i).area=sum([SD.area]);
    elseif i < min(env.cat_Q)
        num=find([attD.Region3]==i); % which shapes to pull
        geom(i).area=sum([SD(num).area]); 
    else % categories, not regions
        num=find([attD.Category]==i-(min(env.cat_Q)-1)); % which shapes to pull
        geom(i).area=sum([SD(num).area]); 
    end
    geom(i).abbrev=labels{i}; % double check
end

%% loop for water fractions -preprocess
% S=union([S.shp]);
% for i= 1:length(F)%Q
%    F(i).shp=intersect(polyshape(F(i).X, F(i).Y), S);
% end
% F=intersect([F.shp], [S.shp]); % <----------------HERE: crop/intersect to study domain
% for i= 1:length(F)%Q
%    F(i).shp=polyshape(F(i).X, F(i).Y);
% end
%% loop for water fractions 

for i=Q
    if i==1
        num_water=find([attWCC.Region3]>0); % which shapes to pull
        num_lakes=find([attF.Region3]>0); % which shapes to pull
    elseif i < min(env.cat_Q)
        num_water=find([attWCC.Region3]==i); % which shapes to pull
        num_lakes=find([attF.Region3]==i); % which shapes to pull
    else % categories, not regions
        num_water=find([attWCC.Category]==i-(min(env.cat_Q)-1)); % which shapes to pull
        num_lakes=find([attF.Category]==i-(min(env.cat_Q)-1)); % which shapes to pull
    end
    geom(i).area_water=sum([attWCC(num_water).Area]); 
    geom(i).area_lakes=sum([attF(num_lakes).Area]); 
    geom(i).count_water_chk=length(num_water); 
    geom(i).count_lakes_chk=length(num_lakes); 
    geom(i).fraction_water=geom(i).area_water/geom(i).area; 
    geom(i).fraction_lakes=geom(i).area_lakes/geom(i).area;
end


%% save
save([saved_dir,'Geom.mat'], 'geom')
fprintf('Mat file saved to %s\n', saved_dir)