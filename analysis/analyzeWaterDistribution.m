% function total=analyzeWaterDistribution(minSize, maxSize)
% function to load and analyze water distribution in aggregate and by
% region.  Makes plots.  Also creates shapefile

% TODO: SDF, other plots?
%% non-function params
clear
minSize=0;
maxSize=1e12; % 1km2
%% params

global plt env % load colors
saveFigs=0;
saveShp=0;
atUCLA=0;

%% directories
if ~isunix
    if atUCLA
        struct_in='J:\output\analysis\distrib.mat';
    else % Brown
        shp_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\shp\distrib.shp';
        struct_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\distrib.mat';
        figs_out='D:\pic\geomFigsBulk\';
        tbl_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\LakeMorphology.xlsx';
        load(env.lake_databases, 'hl_fused')
%         load 'D:\GoogleDrive\Research\Lake distributions\regionLabels2_abrev.mat'
        load(env.labels_in)
    end
else
    struct_in='/Volumes/Galadriel/output/analysis/distrib.mat';
    figs_out='/Volumes/Galadriel/output/pic/geomFigsBulk/';
end

%% load and reshape input data
 load(struct_in); 
%% apply size filters
x1_0=[hl_fused.Area];
p1_0=[hl_fused.Perimeter];
%% regions

%% plot
close all
    % which plots to draw
prelimPlots=0;
plotRegions=0; % plot all regions or just total
plotArea=1;
plotAreaLine=0;
plotCumArea=1;
plotPerim=0;
plotCumPerim=0;
plotDevel=0;
plotAreaDevel=0;
plotBins=0;
plotFit=0; % plot pareto fit
plotMacDonald=0;
plotVerpoorter=0; % must also turn on plotBins
regionsStatus=3; % 1 for pre-2/26, 2 for finer grained, 3 for "final" (post 3/11)

if plotRegions
    i_end=length(regions);
else
    i_end=1;
end

try
    if regionsStatus==2
        for i=1:length(hl_fused)
            hl_fused(i).Region=hl_fused(i).Region2;
        end
    elseif regionsStatus==3
        for i=1:length(hl_fused)
            hl_fused(i).Region=hl_fused(i).Region3_1;
        end
    end
        hl_fused=rmfield(hl_fused, {'Region_1', 'Region2_1', 'Region3_1'});
catch 
    warning('Error.')
end
% ev=sort(10-logspace(-4, 1, 200)); % edge vector to use for binning
ev_ar=logspace(-4, 2, 200);
ev_per=logspace(log10(0.016), 0, 200);
% ev=1e-4:0.002:10;
close all
for i=[1, env.regions_Q, env.cat_Q] %1:13 %[1:8,12,13] % i is number of regions
    % subset by region % x1_0 was all, x1 is just region
    if i==1
        x1=x1_0;
        p1=p1_0;
    elseif i==21 %yukon flats
        x1=[yf.Area];
        p1=[yf.Perimeter];
    elseif i >= 22 % potholes    etc                                                     75
         x1=x1_0([hl_fused.(env.category)]==(i-(min(env.cat_Q)-1)));
         p1=p1_0([hl_fused.(env.category)]==(i-(min(env.cat_Q)-1)));
    else
        x1=x1_0([hl_fused.(env.region)]==(i));
        p1=p1_0([hl_fused.(env.region)]==(i));
    end
    
        % SDF
    SDF=p1./(2.*sqrt(pi*x1));
    if i==1
        SDF_all=p1./(2.*sqrt(pi*x1));
    end
        % count total land and water and other stats
        % percent under 0.01
    g(i).region=labels{i};    
    g(i).ar=x1;
    g(i).ar01=g(i).ar(g(i).ar<0.01); % 1 ha
    g(i).ar001=g(i).ar(g(i).ar<0.001); %< 0.001 km2 or 0.1 ha
    g(i).ar0001=g(i).ar(g(i).ar<0.0001); 
    
    g(i).per=p1;
    g(i).per01=g(i).per(g(i).ar<0.01); % 1 ha
    g(i).per001=g(i).per(g(i).ar<0.001); %< 0.001 km2 or 0.1 ha
    g(i).per0001=g(i).per(g(i).ar<0.0001); 
    
%     total(i).land=sum([abun(abun_msk).land]);
%     total(i).water=sum([abun(abun_msk).water]);
%     total(i).area=total(i).land + total(i).water;
%     total(i).lim=total(i).water/(total(i).water+total(i).land);
%     total(i).lim2=sum([abun(abun_msk).lim].*([abun(abun_msk).land]+...
%         [abun(abun_msk).water])/sum([abun(abun_msk).land]+...
%         [abun(abun_msk).water])); % double check...
    total(i).region=labels{i};
    total(i).count=length(x1);
    total(i).minSize=minSize;
    total(i).maxSize=maxSize;
    total(i).perUnder001=length(g(i).ar001)/total(i).count;
    total(i).perUnder0001=length(g(i).ar0001)/total(i).count;
    total(i).perUnder01=length(g(i).ar0001)/total(i).count;
            % area percentage
    total(i).ArPerUnder01=sum(g(i).ar01)/sum(g(i).ar); 
    total(i).ArPerUnder001=sum(g(i).ar001)/sum(g(i).ar); 
    total(i).ArPerUnder0001=sum(g(i).ar0001)/sum(g(i).ar);
    
    total(i).PerimPerUnder01=sum(g(i).per01)/sum(g(i).per);
    total(i).PerimPerUnder001=sum(g(i).per001)/sum(g(i).per);
    total(i).PerimPerUnder0001=sum(g(i).per0001)/sum(g(i).per);
    
    total(i).MeanArea=mean([g(i).ar]);
    total(i).MedArea=median([g(i).ar]);
    total(i).MeanPerim=mean([g(i).per]);
    total(i).MedPerim=median([g(i).per]);
    total(i).MeanSDF=mean([SDF]);
    total(i).MedSDF=median([SDF]);
    total(i).AreaSD=std([g(i).ar]);
    total(i).PerSD=std([g(i).per]);
    
    pd(i)=fitdist(x1(:), 'GeneralizedPareto', 'Theta', 0.99*minSize/1e6);
    lnd(i)=fitdist(x1(:), 'Lognormal');
    total(i).a=pd(i).sigma; % size param
    total(i).c=pd(i).k; % shape param
    total(i).k=pd(i).theta;
    total(i).logn_mu=lnd(i).mu; % size param
    total(i).logn_sigma=lnd(i).sigma; % shape param
    
  
    par_cdf{i}=gpcdf(ev_ar,pd(i).k,pd(i).sigma,pd(i).theta);
    par_pdf{i}=gppdf(ev_ar,pd(i).k,pd(i).sigma,pd(i).theta);
%     par_cdf_per{i}=gpcdf(ev_per,pd(i).k,pd(i).sigma,pd(i).theta);
%     par_pdf_per{i}=gppdf(ev_per,pd(i).k,pd(i).sigma,pd(i).theta);    
    % plot histogram of area stats
    %       figure
    %     [N,edges] = histcounts(X,edges)
    
        % area
    if plotArea && i==1
    figure%(1)
        histogram(x1, ev_ar, 'FaceColor','auto', 'Normalization', 'count'); xlabel('Area ($km^2$)'); ylabel('Count');
        title({'Area distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
        set(gca, 'YScale', 'log', 'XScale', 'log')
%         annotation(gcf,'textbox',...
%             [0.72 0.70 0.25 0.16],...
%             'String',['n = ', num2str(sum(rshp_msk(:,i)))],...
%             'LineStyle','none',...
%             'FontSize',19,...
%             'FitBoxToText','off');
        xlim([0 10])
        grid on; box on
        if plotFit
        hold on
            plot(ev_ar, par_pdf{i})
            legend({labels{i}, 'Pareto PDF fit'}, 'location', 'best')
        hold off
        end
    end
    
    if plotCumArea && i==1% only make these plots for total extent  
            % cumulative area
%             ev_ar=logspace(-4, 5, 30);
            figure; hold on
            [N,edg]=histcounts(x1, ev_ar, ...
                'Normalization', 'cumcount');
            h=plot(edg(1:end-1), max(N)-N, 'k');
            set(gca, 'LineWidth', 1)
            xlabel('Area ($km^2$)'); ylabel('Number of lakes greater than given area');
            
                % add vert line
            hold on; plot([0.35 0.35], [1e2 1e5],'-r'); hold off
            axis([0.0001            max(ev_ar)          100        1e+05]);
            
%         figure(2); hold on
            

            
                % plot pareto fit
            if plotFit
                hold off
                [CdfY,CdfX] = ecdf(g(i).ar(rshp_msk(:,i))/1e6,'Function','survivor'); 
                plot(CdfX, CdfY)
                hold on
                YPlot = cdf(pd(i),ev_ar);
                YPlot = 1 - YPlot;
                hLine = plot(ev_ar,YPlot)
                legend({'AirSWOT extent', 'Pareto CDF fit'}, 'location', 'best')
            end
            hold off
            xlabel('Area ($km^2$)'); ylabel('Number of lakes of greater area');
%             title({['Cumulative Area distribution (lakes ', num2str(minSize),' - ',num2str(maxSize),' $km^2$)']})%, ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            grid on; box on
            if i==i_end % last time
                if plotRegions
                    legend(labels, 'location', 'best', 'FontSize', 15)
                else
                end
            end
    end
    if plotAreaLine && i==1% only make these plots for total extent  
            % pdf area
%             ev_ar=logspace(-4, 0, 50);

            figure; hold on
            [N,edg]=histcounts(x1, ev_ar, ...
                'Normalization', 'countdensity');
            box on
            grid on
%             plot(edg(1:end-1), max(N)-N);
            plot(edg(1:end-1), N)
            
            xlabel('Area ($km^2$)'); ylabel('Number of lakes of greater area');
%             title({['Area distribution (lakes ', num2str(minSize),' - ',num2str(maxSize),' $km^2$)']})%, ['region: ', labels{i}]}, 'Interpreter', 'none')
            
            set(gca, 'YScale', 'log', 'XScale', 'log')
    end
    if plotPerim
            % perim
        figure%(3)
            h=histogram(g(i).per(rshp_msk(:,i))/1e3, ev_ar, 'FaceColor','auto'); xlabel('Perimeter (km)'); ylabel('Count');
            title({'Perimeter distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(sum(rshp_msk(:,i)))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
            xlim([0.01 10])
    end
    
        if plotCumPerim % only make these plots for total extent  
            % cumulative area
                    figure(4); hold on
            [N,edg]=histcounts(g(i).per(rshp_msk(:,i))/1e3, ev_per, ...
                'Normalization', 'cumcount');
            plot(edg(1:end-1), max(N)-N); xlabel('Perimeter km)'); ylabel('Count of lakes greater than given perimeter');
            
            
%         figure(2); hold on
            

            
                % plot pareto fit
            if plotFit
                hold off
                [CdfY,CdfX] = ecdf(g(i).per(rshp_msk(:,i))/1e3,'Function','survivor'); 
                plot(CdfX, CdfY)
                hold on
                YPlot = cdf(pd(i),ev_per);
                YPlot = 1 - YPlot;
                hLine = plot(ev_per,YPlot)
                legend({'AirSWOT extent', 'Pareto CDF fit'}, 'location', 'best')
            end
            hold off
            xlabel('Perimeter (km)'); ylabel('Number of lakes of greater perimeter');
            title({['Cumulative Perimeter distribution (lakes ', num2str(minSize),' - ',num2str(maxSize),' $km^2$)']})%, ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            if i==i_end % last time
                if plotRegions
                    legend(labels, 'location', 'best', 'FontSize', 15)
                else
                end
            end
    end
    
    if plotDevel
                % perim/area
        figure%(4)
            h=histogram(abun_rshp.SDF(rshp_msk(:,i)), 'FaceColor','auto'); xlabel('SDI'); ylabel('Count');
            title({'Shoreline Development Index' , ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'linear')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(length(SDF))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
    end
    
    if plotAreaDevel
                % area vs SDI
        figure%(5)
            plot(g(i).ar(rshp_msk(:,i))/1e6, abun_rshp.SDF(rshp_msk(:,i)), '.')
            ylabel('SDI'); xlabel('Area ($km^2$)');
            title({'SDI vs area' , ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(sum(rshp_msk(:,i)))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
    end

end
regions={total(~isnan([total.MeanArea])).region}; % list of labels
if prelimPlots % additional summary plots
    for i=1:length(total)
       if isempty(total(i).MeanArea)
        total(i).MeanArea=NaN;
       end 
    end
%     figure
%     histogram([abun.lim]); title('Distribution of water fractions')
%     xlabel('Limnicity per ABoVE tile (\%)')
%     ylabel('Count')

    
% %     figure
% %     bar([total.lim]*100)
% %     set(gca, 'XTickLabel', {total.region}, 'XTickLabelRotation', 90)
% %     title('Water fraction by region')
    
    %
    
    figure
    b(1)=bar([total.perUnder001]*100, 'FaceColor', 'flat');
    set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
    title('Percent of lakes under 0.001 $km^2$')
    ylabel('Percent')
%     b.CData(2,:) = [0 0.8 0.8];
    
    figure
    b(2)=bar([total.ArPerUnder001]*100, 'FaceColor', 'flat');
    set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
    title('Percent of areas under 0.001 $km^2$')
    ylabel('Percent')
    
    figure
    b(3)=bar([total.PerimPerUnder001]*100, 'FaceColor', 'flat');
    set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
    title('Perimeters from lakes under 0.001 $km^2$')
    ylabel('Percent')
    
%     figure
%     bar([total.perUnder0001]*100)
%     set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
%     title('Percent of lakes under 0.0001 $km^2$')
    
%     figure
%     bar([total.ArPerUnder0001]*100)
%     set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
%     title('Percent of areas under 0.0001 $km^2$')

%     figure
%     bar([total.PerimPerUnder0001]*100)
%     set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
%     title('Percent of perimeters from lakes under 0.0001 $km^2$')
    
%     figure
%     bar([total.PerimPerUnder01]*100)
%     set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
%     title('Percent of perimeters from lakes under 0.01 $km^2$')
    
%     boxArea=[];
%     boxRegion='';
%     for k=1:length(g)
%         boxArea=[boxArea; [g(k).ar]'];
%         boxRegion=[boxRegion; repmat(g(k).region, length([g(k).ar]),1)];
% %         boxArea(k,:)=[g(k).ar];
%     end
% %     boxRegion=cellstr(boxRegion);
%     figure
%     b(4)=boxplot(boxArea, boxRegion);
    
    
    figure
    b(4)=bar([total.MeanArea], 'FaceColor', 'flat');
    set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
    title('Median lake area')
    ylabel('($km^2$)')
%     hold on
%     errhigh=[total.MedArea]+[total.AreaSD];
%     errlow=[total.MedArea]-[total.AreaSD];
%     x=1:length([total.MedArea]);
%     data=[total.MedArea];
%     er = errorbar(x,data,errlow,errhigh);    
%     er.Color = [0 0 0];                            
%     er.LineStyle = 'none';  
%     hold off
    
    
    figure
    b(5)=bar([total.MedPerim], 'FaceColor', 'flat');
    set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
    title('Median lake perimeter')
    ylabel('($km$)')
    
    figure;
    b(6)=scatter(g(1).ar, SDF_all);
    xlabel('Area ($km^2$)'); ylabel('SDF')
    set(gca, 'YScale', 'log', 'XScale', 'log')
    
    figure
    b(7)=bar([total.MedSDF], 'FaceColor', 'flat');
    set(gca, 'XTickLabel', regions, 'XTickLabelRotation', 90)
    title('Median SDF')
    
    figure
    b(8)=scatter([total.MedArea], [total.MedSDF]);
    xlabel('Area ($km^2$)'); ylabel('SDF')
    box on
    

if plotMacDonald %( 3 is Yukon flats)
    i=3;
    % here
    ev_mac=logspace(1.5, 7, 40);
    histogram(g(i).ar(rshp_msk(:,i)), ev_mac, 'FaceColor','auto', 'Normalization', 'count'); xlabel('Area ($m^2$)'); ylabel('Count');
        title({'Area distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
        set(gca, 'YScale', 'lin', 'XScale', 'log')
        annotation(gcf,'textbox',...
            [0.72 0.70 0.25 0.16],...
            'String',['n = ', num2str(sum(rshp_msk(:,i)))],...
            'LineStyle','none',...
            'FontSize',19,...
            'FitBoxToText','off');
        xlim([0 1e7])
        if plotFit
        hold on
            plot(ev_ar, par_pdf{i})
            legend({labels{i}, 'Pareto PDF fit'}, 'location', 'best')
        hold off
        end
        

end

%% bin perimeter by area

    
if plotBins 
    fn=get(gcf,'Number');
    plotLogspace=1;
    plotSemiLogSpace=0;
        % params
    if plotLogspace
%         edges1_area=[4e-5 1e-4 1e-3 0.01 0.1 1];
%         edges1_area=[ min(g(i).ar) 4e-3 0.04 0.4 max(g(i).ar)+1];
        edges1_area=logspace(log10(min(x1_0)), log10(max(x1_0)), 50);
%         edges1_perim=[min(g(i).per) 0.19 1.9 19  max(g(i).per)+1];
        edges1_perim=logspace(log10(min(p1_0)), log10(max(p1_0)), 50);        
    else
        if 1==0 %plotVerpoorter
%             edges1_area=logspace(-4, 1, 6);
%             edges1_perim=
        else
            edges1_area=linspace((min(x1_0)), log10(max(x1_0)), 50);
            edges1_perim=linspace((min(p1_0)), log10(max(p1_0)), 50);
        end
    end
    
    
        % just area
    
%     figure
%     subplot(211)
    [counts, ~, binsArea] = histcounts(x1_0,edges1_area);
%     histogram(g(i).ar,edges1_area);
%     if plotLogspace
%         set(gca, 'YScale', 'lin', 'XScale', 'log')
%     elseif plotSemiLogSpace
%         set(gca, 'YScale', 'log', 'XScale', 'lin')
%     end
%     xlabel('Area ($km^2$)'); ylabel('count')
%     title('Area histogram')
    
%     subplot(212)
    aggAreabyArea=zeros(size(counts)); % init
    for i=1:max(binsArea)
        aggAreabyArea(i)=sum(x1_0(binsArea==i)); 
    end
    if sum(binsArea==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsArea==0))
    end
    figure
    b(fn+1)=histogram('BinEdges',edges1_area, 'BinCounts', aggAreabyArea)
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    elseif plotSemiLogSpace
        set(gca, 'YScale', 'log', 'XScale', 'lin')
    end
    xlabel('Area ($km^2$)'); ylabel('Sum of binned areas ($km^2$)')
    title('Binned by area') 
    
        % just perim
%     figure
%     subplot(211)
%     histogram(g(i).per,edges1_perim);
%     if plotLogspace
%         set(gca, 'YScale', 'lin', 'XScale', 'log')
%     elseif plotSemiLogSpace
%         set(gca, 'YScale', 'log', 'XScale', 'lin')
%     end
%     xlabel('Perimeter ($km$)'); ylabel('count')
%     title('Perimeter histogram')
%         
%     subplot(212)
    figure
    [counts, ~, binsPerim] = histcounts(p1_0,edges1_perim);
    aggPerimbyPerim=zeros(size(counts)); % init
    for i=1:max(binsPerim)
        aggPerimbyPerim(i)=sum(p1_0(binsPerim==i));
        
    end
    if sum(binsPerim==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsPerim==0))
    end
    b(fn+2)=histogram('BinEdges',edges1_perim, 'BinCounts', aggPerimbyPerim)
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    elseif plotSemiLogSpace
        set(gca, 'YScale', 'log', 'XScale', 'lin')
    end
    xlabel('Perim ($km^2$)'); ylabel('Sum of binned perimeters (km)')
    title('Binned by perimeter')
    
        % bin perim by area
    figure    
    aggPerimbyArea=zeros(size(counts)); % init
    for i=1:max(binsArea)
        aggPerimbyArea(i)=sum(p1_0(binsArea==i));
    end
    if sum(binsArea==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsPerim==0))
        end
    b(fn+3)=histogram('BinEdges',edges1_area, 'BinCounts', aggPerimbyArea)
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    elseif plotSemiLogSpace
        set(gca, 'YScale', 'log', 'XScale', 'lin')
    end
    xlabel('Area ($km^2$)'); ylabel('Sum of binned perimeters (km)')
    title('Binned by area')
end

%% verpoorter plots
% close all
plotLogspace=0;
plotSideways=0;
if plotVerpoorter
    if plotLogspace
%         edges1_area=logspace(log10(5e-5), log10(5), 6);
        edges1_area=logspace(-4, 0, 5);
    else
        edges1_area=linspace(5e-5, 1, 50);
    end
    [counts, ~, binsArea] = histcounts(x1_0,edges1_area);    
    aggPerimbyArea=zeros(size(counts)); % init
    aggAreabyArea=zeros(size(counts)); % init
    aggCountbyArea=zeros(size(counts)); % init
    for i=1:max(binsArea)%:-1:1
        aggPerimbyArea(i)=sum(p1_0(binsArea==i))/1000;
        aggAreabyArea(i)=sum(x1_0(binsArea==i))/1e6;
        aggCountbyArea(i)=length(x1_0(binsArea==i));
    end
    if sum(binsArea==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsPerim==0))
    end
    
    figure
    histogram('BinEdges',edges1_area, 'BinCounts', aggPerimbyArea)
    if plotLogspace==1
        set(gca, 'YScale', 'lin', 'XScale', 'log')
        set(gca, 'XTick', edges1_area,'view',[90 -90])
    elseif plotLogspace==0
        set(gca, 'YScale', 'lin', 'XScale', 'lin')
    end
    xlabel('Area ($km^2$)'); ylabel('Sum of binned perimeters (km)')
    title('Perimeter histogram binned by area')
    if plotSideways
    set(gca,'view',[90 -90])
    end
    figure
    histogram('BinEdges',edges1_area, 'BinCounts', aggAreabyArea)
    if plotLogspace==1
        set(gca, 'YScale', 'lin', 'XScale', 'log')
        set(gca, 'XTick', edges1_area,'view',[90 -90])
    elseif plotLogspace==0
        set(gca, 'YScale', 'lin', 'XScale', 'lin')
    end
    xlabel('Area ($km^2$)'); ylabel('Sum of binned areas (km2)')
    title('Area histogram binned by area')
    if plotSideways
    set(gca,'view',[90 -90])
    end
    
    figure
    histogram('BinEdges',edges1_area, 'BinCounts', aggCountbyArea)
    if plotLogspace==1
        set(gca, 'YScale', 'log', 'XScale', 'log')
        set(gca, 'XTick', edges1_area,'view',[90 -90])
    elseif plotLogspace==0
        set(gca, 'YScale', 'lin', 'XScale', 'lin')
    end
    xlabel('Area ($km^2$)'); ylabel('Number of lakes within bin')
    title('Histogram binned by area')
    if plotSideways
    set(gca,'view',[90 -90])
    end
end

%% touch up figs
    for d=1:get(gcf, 'Number')
       figure(d)
       box on
       set(gca, 'LineWidth', 1)
       if d~=6 && d~=8 && d<9
           b(d).CData=repmat(plt.c, [size(b(d).CData, 1),1, 1]);
           b(d).CData(1,:)=[0.71 1 1];
       else
           try
                b(d).CData=plt.c;
           catch
               b(d).FaceColor=plt.c;
           end
       end
    end
end

%% output shapefile

if saveShp
    disp('Saving shape...')
    S=mappoint(abun_rshp.long,abun_rshp.lat, abun_rshp);
    shapewrite(S, shp_out);
end


%% computations

% load('D:\GoogleDrive\Research\Lake distributions\LakeDatabases.mat', 'hl_global');
% lim=0.5; %0.34
% hla=[hl.Lake_area];
% sum(hla(hla<lim))/sum(hla)*100 %5.5

% load('D:\GoogleDrive\Research\Lake distributions\LakeDatabases.mat', 'hl_global')
% g(1).ar35=g(1).ar(g(1).ar<0.35); 
% g(1).ArPerUnder35=sum(g(1).ar35)/sum(g(1).ar);
% g(1).ar_glob=[hl_global.Lake_area];
% g(1).ar35_glob=g(1).ar_glob(g(1).ar_glob<0.35); 
% g(1).ArPerUnder35_glob=sum(g(1).ar35_glob)/sum(g(1).ar_glob);

%% output stats table

% tbl=struct2table(total);
% writetable(tbl, tbl_out);  % TODO: load regional fitting data

%% save mat file
% load('D:\GoogleDrive\Research\Lake distributions\savedData\tmp\analyzeWaterDistribution_old_temp.mat', 'total_lim')
for i=1:length(total)
    total(i).lim=[];
end
save(env.analyzeWaterDistribution, 'g', 'total', 'regions', 'SDF_all')