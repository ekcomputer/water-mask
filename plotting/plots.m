% for plotting illustrative examples.  Uses namespace of after running
% BP_OBIA_Devel.m

pic_dir='D:\pic\AGU2018Figs\';
zoom_bounds= [2036.8       2511.2       1417.4       1813.4];

figure
imagesc(cir); axis image

%% loop
f=[13, 3, 11, 10, 12] % fig numbers
descript={'Image', 'Mean', 'MeanSP','TextureSP', 'Final'};
for i=1:length(f)
    figure(f(i))
    axis image
    set(gca, 'Visible', 'off')
    set(gca, 'Position', [0 0 1 1])
    set(gcf,'Resize','off')
    set(gcf, 'InnerPosition', [-1487.5 -532.5 1227  1023], 'OuterPosition',...
        [-1487.5 -532.5 1107  1023])
%     set(gca, 'OuterPosition', [0, 0, 3176/2651, 1])
%     set(gca, 'Position', [0, 0, 3176/2651, 1])
%     set(gcf, 'OuterPosition', get(gcf, 'Position'))
%     set(gca, 'Position', [0 0 1.1 1.1])
    name_out=[pic_dir, name_in(1:end-4), '_', descript{i}, '.jpg'] 
        % save
    saveas(gcf, name_out)
    axis(zoom_bounds)
    name_out=[pic_dir, name_in(1:end-4), '_', descript{i}, '_zoom', '.jpg'] 
    
        %save
    saveas(gcf, name_out)
end


%%
% print(gcf,name_out,'-djpeg')
% saveas(gcf, name_out)