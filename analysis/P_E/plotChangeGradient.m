% inputs text file output from ArcGIS with lake area change stats, then
% plots based on lat/long

clear; close all

f.tbl_in='J:\output\lakeChangeStats\UpperAthaChange.txt';
chg=table2struct(readtable(f.tbl_in));
chg=chg([[chg.pchange]<9000]);

    % fig 1
plot([chg.lat], [chg.pchange], '.', 'MarkerSize', 12)

axis([54.1 55.7 -500 2000]);
p=polyfit([chg.lat], [chg.pchange], 1)
r2=corrcoef([chg.lat], [chg.pchange]).^2
hold on; fplot(@(x) p(1)*x + p(2)); hold off
ylabel('Area Change (\%)')
xlabel('Latitude')
annotation(gcf,'textbox',...
    [0.63 0.72 0.25 0.16],...
    'String',{['r^2 = ', num2str(r2(3))],'slope = ',num2str(p(1))},...
    'LineStyle','none',...
    'FontSize',15,...
    'FitBoxToText','off');

    % fig 2
figure;
plot([chg.lat], [chg.change], '.', 'MarkerSize', 12)

axis([54.1 55.7 -5000 50000]);
p=polyfit([chg.lat], [chg.change], 1)
r2=corrcoef([chg.lat], [chg.change]).^2
hold on; fplot(@(x) p(1)*x + p(2)); hold off
% ylim([-800 800]);
ylabel('Area Change ($m^2$)')
xlabel('Latitude')
annotation(gcf,'textbox',...
    [0.63 0.72 0.25 0.16],...
    'String',{['r^2 = ', num2str(r2(3))],'slope = ',num2str(p(1))},...
    'LineStyle','none',...
    'FontSize',15,...
    'FitBoxToText','off');