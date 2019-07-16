% uses new toolbox
addpath D:\Dropbox\Matlab\DownloadedCode\PowerLawAaronClauset\powerlaws_full_v0.0.10-2012-01-17\powerlaws_full_v0.0.10-2012-01-17
clear

%% params

prePlots=false;
multiPlots=true;
singlePlots=false;
plotPDFs=true;
plotCDFs=false;
old=false; % updated data
%% load data

if old
    compareDistribs
    % load vars
    le=length(abun_rshp.ar);
    % msk1=le-unique(sort(round(logspace(0, 4.873, 500)), 'descend'));
    x1=abun_rshp.ar(msk1);
    x2=glwd.ar;
else
    load('D:\GoogleDrive\Research\Lake distributions\LakeDatabases.mat', 'hl_fused','dcs','hl')
%     load('D:\GoogleDrive\Research\Lake distributions\AlaskaLakes\Yukonflats.mat')
    le=length([dcs.Area]);
    msk1=unique(randi(le, [1000,1]));
%     x1_0=[dcs.Area]/1e6; %x1=x1(msk1);
    x1_0=[hl_fused.Area]/1e1; 
    x2=[hl.Lake_area];
    x3=[dcs.Area];
end
load 'D:\GoogleDrive\Research\Lake distributions\regionLabels.mat'




%% fit power laws
% [alpha1, xmin1, L1]=plfit(x1_0)
% [alpha2, xmin2, L2]=plfit(x2)
% [alpha3, xmin3, L3]=plfit(x3)
    % save time by using pre-run values:
alpha1=1.7867; xmin1=0.0076324; L1 = 5589.5;
alpha2=1.84; xmin2=1.7; L2 = -705.42;
alpha3 = 1.4387; xmin3 = 0.000145;L3 =   2.1847e+05;
%% plot power laws - CDF
if prePlots
    plplot(x1_0, xmin1, alpha1); title ('DCS-HydroLAKES-fused')
    plplot(x2, xmin2, alpha2); title('HydroLAKES')
    plplot(x3, xmin3, alpha3); title ('DCS')
end
%% plot power laws - PDF
close all
Q=[1:5, 10:13,22];%[1];% 12 9] %1:13
Q=1;
c=1; % counter
for region=Q
        % subset by region % x1_0 was all, x1 is just region
    if region==1
        x1=x1_0;
    elseif region==21 %yukon flats
        x1=[yf.Area];
    elseif region==22 % shield w and e
        x1=x1_0([hl_fused.Region]==7 | [hl_fused.Region]==8);
    else
        x1=x1_0([hl_fused.Region]==region);
    end
        % compute fits for region
    if region==1
        [alpha(region), xmin(region), L(region)]=deal(alpha1, xmin1, L1);
    else
        [alpha(region), xmin(region), L(region)]=plfit(x1);
    end
        % plot power law CDF
    if plotCDFs && singlePlots
        plplot(x1, xmin(region), alpha(region)); title({'Area distribution (fused)', labels{region}}, 'Interpreter', 'none')
        annotation(gcf,'textbox',...
        [0.68 0.70 0.25 0.16],...
        'String',{['\alpha = ', num2str(alpha(region))],['n = ', num2str(length(x1))],...
        ['x\_min = ', num2str(xmin(region))]},'LineStyle','none','FontSize',19,'FitBoxToText','off');
    end
        % plot power law PDF
    if plotPDFs && singlePlots
        ev_ar=logspace(-5, 4, 100);
        [x1_counts, x1_bins]=histcounts(x1, ev_ar, 'Normalization', 'countdensity');
        x1_bins=movmean(x1_bins, 2, 'Endpoints', 'discard');
        figure; plot(x1_bins, x1_counts, 'bo', 'MarkerSize', 8)
        title({'Area distribution', labels{region}}, 'Interpreter', 'none')
        xlabel('Area ($km^2$)'); ylabel({'Normalized number of lakes', 'within size class'});
        set(gca, 'YScale', 'log', 'XScale', 'log')
        annotation(gcf,'textbox',...
        [0.72 0.70 0.25 0.16],...
        'String',['n = ', num2str(length(x1))],...
        'LineStyle','none','FontSize',19,'FitBoxToText','off');
        xlim([0 100])
        set(gca, 'YScale', 'log', 'XScale', 'log')
    end
    if multiPlots
        Markers = {'+','o','*','x','v','diamond','^','square','>','<'};
        if plotCDFs
            figure(3); hold on
            plplot_simple(x1, xmin(region), alpha(region));
            hold off
            if region==Q(end)
                [h, ~, plots] =legend(labels(Q), 'location', 'best', 'interpreter', 'tex')
                    % set legend text to match line color
                for idx = 1:length(h.String)
                    h.String{idx} = ['\color[rgb]{' num2str(plots(idx).Color) '} ' h.String{idx}]
                end  

                    % set different markers
                c=get(gca, 'Children');
                for j=1:length(c)
                    set(c(j), 'Marker', Markers{j}, 'MarkerSize', 2);
                end
            end
        end
        
        % multi PDFs
        
        if plotPDFs
            figure(4); hold on
            ev_ar=logspace(-5, 4, 100);
            [x1_counts, x1_bins]=histcounts(x1, ev_ar, 'Normalization', 'countdensity');
            x1_bins=movmean(x1_bins, 2, 'Endpoints', 'discard');
            plot(x1_bins, x1_counts, Markers{c}, 'MarkerSize', 2); drawnow
            title({'Area distribution', labels{region}}, 'Interpreter', 'none')
            xlabel('Area ($km^2$)'); ylabel({'Normalized number of lakes', 'within size class'});
            set(gca, 'YScale', 'log', 'XScale', 'log')
            xlim([0 100])
            set(gca, 'YScale', 'log', 'XScale', 'log')
            hold off
            if region==Q(end)
                [h, ~, plots] =legend(labels(Q), 'location', 'best', 'interpreter', 'tex')
                        % set legend text to match line color
                    for idx = 1:length(h.String)
                        h.String{idx} = ['\color[rgb]{' num2str(plots(idx).Color) '} ' h.String{idx}]
                    end  
            end
        end
    end
    c=c+1;
end
%% goodness of fits
% [p1, gof1] = plpva(x1, xmin1)  % p=0.5?
% [p2, gof2] = plpva(x2, xmin2)