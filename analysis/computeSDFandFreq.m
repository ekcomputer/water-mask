% computes shoreline development index (SDF) for stats list
% also computes frequency in count/100 km28

% run this after waterDistribution and appendWaterDistribution.m
clear
struct_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\distrib.mat';
load(struct_in);
for i=1:length(abun)
    for j=1:length(abun(i).stats)
                abun(i).stats(j).SDF=abun(i).stats(j).Perimeter/(2*sqrt(pi*abun(i).stats(j).Area));
    end
    abun(i).freq_min40=length(abun(i).stats)/(abun(i).water+abun(i).land)*100e6;
end

% save(struct_in, 'abun')