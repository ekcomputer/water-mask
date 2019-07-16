% script to plot global lake area estimates as bar chart

y=[2.803
3.2
4.46
3.8
3
5
4.95
]
str={
    'Maybeck 1995'
'Lehner and Doll 2004'
'Downing et al. 2006'
'Downing et al. 2012'
'MacDonald et al. 2012'
'Raymond et al. 2013'
'Verpoorter et al 2014'
'Holgerson and Raymond 2016'
'Allen and Pavelsky 2018'

}


hB=bar(y);          % use a meaningful variable for a handle array...
hAx=gca;            % get a variable for the current axes handle
hAx.XTickLabel=str; % label the ticks
hT=[];              % placeholder for text object handles
for i=1:length(hB)  % iterate over number of bar objects
  hT=[hT text(hB(i).XData+hB(i).XOffset,hB(i).YData,num2str(hB(i).YData.','%.1f'), ...
                          'VerticalAlignment','bottom','horizontalalign','center',...
                          'FontSize', 19)];
end
ylabel('Area (million $km^2$)')
set(gca, 'XTickLabelRotation', 40)
title('Global Lake Area Estimates', 'FontSize', 38)
ylim([0 5.5])