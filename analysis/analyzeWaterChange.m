% script to plot water area change stats
% option to recalculate percent change (divide by max, not first)

clear

%% params
reprocess=0; % don't reprecoess and write new file
minSize=4000;
%% dir and load
pth_in='J:\Final\analysis\areaChange\shp\aggregate\water_change.shp';
pth_out='J:\Final\analysis\areaChange\shp\aggregate\water_change_2.shp';
shp=mappoint(shaperead(pth_in));

 %% reprocess
if reprocess
    shp.pchange3=100*(shp.size_b./max(shp.size_a, shp.size_b)-shp.size_a./max(shp.size_a, shp.size_b));

    % write file
    shapewrite(shp, pth_out)
end

%% pre-process

pchange=shp.pchange3(shp.pchange3<100 & shp.pchange3>-100);
pchange_mskd=shp.pchange3(shp.size_a>minSize | shp.size_b>minSize);
%% plot stats

figure(1)
histogram(pchange, [-100:100])
% set(gca, 'XScale', 'log')
xlabel('\% change'); ylabel('Count'); title('Lake area change')
    annotation(gcf,'textbox',...
        [0.73 0.75 0.25 0.16],...
        'String',['n = ', num2str(length(pchange))],...
        'LineStyle','none',...
        'FontSize',19,...
        'FitBoxToText','off');
    