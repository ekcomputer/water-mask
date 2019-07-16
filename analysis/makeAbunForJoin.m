% preps data to join in ArcGIS
% runs analyzeWaterDistribution.m first
analyzeWaterDistribution
% vars
lookup_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\join.xlsx';
fact=100;% set to 1 to view percents as fractions
% loop
for i=1:length(abun)
    abun(i).perUnder01=sum([abun(i).stats.Area]<0.01e6)/length(abun(i).stats)*fact;
    abun(i).perUnder001=sum([abun(i).stats.Area]<0.001e6)/length(abun(i).stats)*fact;
    abun(i).perUnder0001=sum([abun(i).stats.Area]<0.0001e6)/length(abun(i).stats)*fact;
    
    ars=[abun(i).stats.Area];
    abun(i).ArPerUnder01=sum(ars<0.01e6)/sum(ars)*fact;
    abun(i).ArPerUnder001=sum(ars<0.001e6)/sum(ars)*fact;
    abun(i).arPerUnder0001=sum(ars<0.0001e6)/sum(ars)*fact;
    
    pers=[abun(i).stats.Perimeter];
    abun(i).PerPerUnder01=sum(pers(ars<0.01e6))/sum(pers)*fact;
    abun(i).PerPerUnder001=sum(pers(ars<0.001e6))/sum(pers)*fact;
    abun(i).PerPerUnder0001=sum(pers(ars<0.0001e6))/sum(pers)*fact;
    
    % change to arc index
    abun(i).file_idx=abun(i).file_idx-1;
    
    % change lim
    abun(i).lim=abun(i).lim*fact;
end

%% save
lookup=rmfield(abun, 'stats');
writetable(struct2table(lookup), lookup_out)
