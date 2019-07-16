% runs analyzeWaterdistrib to test min size param
% plots pareto params vs min size

clear
savePlots=0;
saveMat=1;
startNew=1;
%% dir
mat_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\total.mat';
figs_out='D:\pic\geomFigsBulk\';
minSize=[40 100 200 500 1000 5000 10000];
maxSize=[1e5 1e6 1e7];

%% run for multiple min/max sizes!
if startNew % recompute
    d=1; % counter
    for i=1:length(minSize)
        for j=1:length(maxSize)
            disp(d)
            total{i,j}.total=analyzeWaterDistribution(minSize(i), 1e6);
             total{i,j}.minSize=minSize(i);
             total{i,j}.maxSize=maxSize(j)
            d=d+1;
        end
    end
else % load from previous run
   load(mat_out) 
end
%% save
if saveMat
save(mat_out, 'total')
end
%% analyze 
    % NOTE: only looks at aggregate stats for whole dataset

for i=1:size(total, 1)
    min(i)= total{i,3}.minSize/1e6;
    c(i)= total{i,3}.total(1).c;
    m(i)= total{i,3}.total(1).logn_mu;
    s(i)= total{i,3}.total(1).logn_sigma;
end
figure
plot(min, c); xlabel('Minimum size ($km^2$)'); ylabel('Pareto c')
title('Fixed max size: 1 $km^2$')

figure
AX=plotyy(min, m, min, s); xlabel('Minimum size ($km^2$)'); ylabel('Lognormal mu')
ylabel(AX(2),'Lognormal sigma')
title('Fixed max size: 1 $km^2$')


clear c
for i=1:size(total, 2)
    max(i)= total{2, i}.maxSize/1e6;
    c(i)= total{2,i}.total(1).c;
end
figure
plot(max, c); xlabel('Maximum size ($km^2$)'); ylabel('Pareto c')
title('Fixed min size: 0.0001 $km^2$')

%% analyze by region
close all
for j=1:length(total{1,1}.total) % for each region
    clear min c
    for i=1:size(total, 1)
        min(i)= total{i,3}.minSize/1e6;
        c(i)= total{i,3}.total(j).c;
    end
    figure(1)
    hold on
    plot(min, c); xlabel('Minimum size ($km^2$)'); ylabel('Pareto c')
%     title({'Fixed max size: 1 $km^2$', ['Region: ',total{1,1}.total(j).region]})
    title({'Fixed max size: 1 $km^2$'})
    if j==length(total{1,1}.total) % last time
            legend({total{1,1}.total.region}, 'location', 'best', 'FontSize', 23)
            set(gca, 'YScale', 'lin', 'XScale', 'log')
    end
    hold off
%     clear min c
%     for i=1:size(total, 2)
%         max(i)= total{2, i}.maxSize/1e6;
%         c(i)= total{2, i}.total(j).c;
%     end
%     figure
%     plot(max, c); xlabel('Maximum size ($km^2$)'); ylabel('Pareto c')
%     title({'Fixed min size: 0.0001 $km^2$', ['Region: ',total{1,1}.total(j).region]})
end

%% save plots
if savePlots
for i=1:get(gcf, 'Number')
    saveas(i, [figs_out, 'ParetoShape_', num2str(i), '.png'])
end

    % save text file pointing to current directory (for this script)
fid=fopen([figs_out, 'Source.txt'], 'w+');
fprintf(fid, '%s\n', pwd);
end