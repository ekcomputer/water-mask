% script to compute areas of each region and water fractions
clear
global env
shp_in=env.shp_in;
load(env.labels_in);
saved_dir=env.saved_dir;
[S,att]=shaperead(shp_in);
Q=[1,env.regions_Q, env.cat_Q];

%% loop
for i= 1:length(S)%Q
   S(i).shp=polyshape(S(i).X, S(i).Y) 
   S(i).area=area(S(i).shp)/1e6;
end

for i=Q
    if i==1
        geom(i).area=sum([S.area]);
    elseif i < 22
        num=find([att.Region3]==i); % which shapes to pull
        geom(i).area=sum([S(num).area]); 
    else % categories, not regions
        num=find([att.Category]==i-21); % which shapes to pull
        geom(i).area=sum([S(num).area]); 
        geom(i).abbrev=labels{i}; % double check
    end
end

%% save
% save([saved_dir,'Geom.mat'], 'geom')