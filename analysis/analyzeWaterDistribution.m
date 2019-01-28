% function total=analyzeWaterDistribution(minSize, maxSize)
% function to load and analyze water distribution in aggregate and by
% region.  Makes plots.  Also creates shapefile

%% non-function params
clear
minSize=50;
maxSize=1e6; % 1km2
%% params

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
    end
else
    struct_in='/Volumes/Galadriel/output/analysis/distrib.mat';
    figs_out='/Volumes/Galadriel/output/pic/geomFigsBulk/';
end

%% load and reshape input data
load(struct_in);
abun_rshp.ar=[]; abun_rshp.per=[]; abun_rshp.lat=[]; abun_rshp.file_idx=[];
abun_rshp.SDF=[]; abun_rshp.LehnerDevel=[]; abun_rshp.long=[];
for i=1:length(abun)
    ar_temp=[abun(i).stats.Area];
    per_temp=[abun(i).stats.Perimeter];
    lat_temp=[abun(i).stats.lat];
    SDF_temp=[abun(i).stats.SDF];
    long_temp=[abun(i).stats.long];
    LehnerDevel_temp=[abun(i).stats.LehnerDevel];
    file_idx_temp=repmat([abun(i).file_idx], [1,length(abun(i).stats)]);
    abun_rshp.ar=[abun_rshp.ar,ar_temp];
    abun_rshp.per=[abun_rshp.per,per_temp];
    abun_rshp.lat=[abun_rshp.lat,lat_temp];
    abun_rshp.file_idx=[abun_rshp.file_idx,file_idx_temp]; % record index for masking
    abun_rshp.SDF=[abun_rshp.SDF,SDF_temp];
    abun_rshp.LehnerDevel=[abun_rshp.LehnerDevel,LehnerDevel_temp];
    abun_rshp.long=[abun_rshp.long,long_temp];
end


%% apply size filters
sizeMsk= ([abun_rshp.ar]>=minSize & [abun_rshp.ar] < maxSize  );
flds=fieldnames(abun_rshp);
for j=1:length(flds)
    abun_rshp.(flds{j})=abun_rshp.(flds{j})(sizeMsk);
end
%% regions

% all
regions{1}=[1:330, -99];
labels{1}='AirSWOT extent';

% north slope
regions{2}=269:272;
labels{2}='North Slope';

% yukon flats
regions{3}=[1+[65	66	67	68	69	70	71	72	73	74	90	91	263	264	265	266	267	268	273	274	275	276	277	278	279	280], -99];
labels{3}='Yukon Flats';

% old crow flats
regions{4}=1+[80	81	82	83	84	85	86	87	88	89];
labels{4}='Old Crow Flats';

% inuvik
regions{5}=1+[55	56	57	58	59	60	61	62	63	64	75	76	77	78	79	92	93	94	95	96	97	98	99	100	101	102	103	104	105	259	260	261	262];
labels{5}='Mackenzie Delta';

% mackenzie river
regions{6}=1+[106	107	108	109	110	111	112	113	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	134	135	136 137 138 139 140 141 142];
labels{6}='Mackenzie River';

% yellowknife W
regions{7}=1+[143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165];
labels{7}='Yellowknife W';

% yellowknife E
regions{8}=1+[44	45	46	47	48	49	50	51	52	53	54	166	167	168	169	170	171	172	173	174	175	176	177	178	179	180	181	182	183	184	185	186	187	188	211	212	213	214	215	216	217	218	219	220	221	222	223];
labels{8}='Yellowknife E';

regions{9}=1+[189:200];
labels{9}='Slave River';

regions{10}=1+[281:287];
labels{10}='Peace-Athabasca Delta';

regions{11}=1+[0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	201	202	203	204	205	206	207	208	209	210	288	289	290	291	292	293	294	295	296	297	298	299	300	301	302	303];
labels{11}='Athabasca River - Edmonton';

regions{12}=1+[26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	224	225	226	227	228	229	230	231	232	233	234	235	236	237	238	239	240	241	242	243	244	245	246	247	248	249	250	251	252	253	254	255	256	257	258	304	305	306	307	308	309	310	311	312	313	314	315	316	317	318	319	320	321];
labels{12}='Edmonton-Saskatoon';

regions{13}=1+[322:329];
labels{13}='North Dakota Pothole Lakes';

%% plot
close all
    % which plots to draw
prelimPlots=0;
plotRegions=1; % plot all regions or just total
plotArea=0;
plotCumArea=0;
plotPerim=0;
plotCumPerim=0;
plotDevel=0;
plotAreaDevel=0;
plotBins=0;
plotFit=0; % plot pareto fit
plotMacDonald=0;

if plotRegions
    i_end=length(regions);
else
    i_end=1;
end
% ev=sort(10-logspace(-4, 1, 200)); % edge vector to use for binning
ev_ar=logspace(-4, 0, 200);
ev_per=logspace(log10(0.016), 0, 200);
% ev=1e-4:0.002:10;
for i=1:i_end % i is number of regions
    regions{i}=intersect(regions{i}, [abun.file_idx]); % use only idx with observations included in complete map
    rshp_msk(:,i)=false(size(abun_rshp.file_idx)); % init; mask selects only files in region of interest
    abun_msk=false(size(abun));
    for k=1:length(abun_msk)
       if find(abun(k).file_idx==regions{i})>0 
           abun_msk(k)=1;
       end
    end
    for j=1:length(regions{i}) % j is number of tiles in each region
        rshp_msk(:,i)=rshp_msk(:,i) | transpose(abun_rshp.file_idx==regions{i}(j)); % true values are regions (files) of interet
    end
    
        % count total land and water and other stats
        % percent under 0.01
    g.ar=[abun_rshp.ar(rshp_msk(:,i))];
    g.ar01=g.ar(g.ar<10000); % 1 ha
    g.ar001=g.ar(g.ar<1000); %< 0.001 km2 or 0.1 ha
    g.ar0001=g.ar(g.ar<100); 
    
    g.per=[abun_rshp.ar(rshp_msk(:,i))];
    g.per01=g.ar(g.per<10000); % 1 ha
    g.per001=g.ar(g.per<1000); %< 0.001 km2 or 0.1 ha
    g.per0001=g.ar(g.per<100); 
    
    total(i).land=sum([abun(abun_msk).land])/1e6;
    total(i).water=sum([abun(abun_msk).water])/1e6;
    total(i).area=total(i).land + total(i).water;
    total(i).lim=total(i).water/(total(i).water+total(i).land);
    total(i).lim2=mean([abun(abun_msk).lim].*([abun(abun_msk).land]+...
        [abun(abun_msk).water])/sum([abun(abun_msk).land]+...
        [abun(abun_msk).water])); % double check...
    total(i).region=labels{i};
    total(i).count=length([abun_rshp.ar(rshp_msk(:,i))]);
    total(i).minSize=minSize;
    total(i).maxSize=maxSize;
    total(i).perUnder001=sum([abun_rshp.ar(rshp_msk(:,i))]<1000)/total(i).count;
    total(i).perUnder0001=sum([abun_rshp.ar(rshp_msk(:,i))]<100)/total(i).count;
    total(i).perUnder01=sum([abun_rshp.ar(rshp_msk(:,i))]<10000)/total(i).count;
            % this takes sum of areas under certain amount and divides by
            % lake:totalwater ratio...
%     total(i).ArPerUnder01=sum(g.ar)/1e6/total(1).water * sum(g.ar01)/1e6; 
%     total(i).ArPerUnder001=sum(g.ar)/1e6/total(1).water * sum(g.ar001)/1e6; 
%     total(i).ArPerUnder0001=sum(g.ar)/1e6/total(1).water * sum(g.ar0001)/1e6; 

    total(i).ArPerUnder01=1/total(1).water * sum(g.ar01)/1e6; 
    total(i).ArPerUnder001=1/total(1).water * sum(g.ar001)/1e6; 
    total(i).ArPerUnder0001=1/total(1).water * sum(g.ar0001)/1e6; 
    
    total(i).PerimPerUnder01=sum(g.per01)/sum(g.per);
    total(i).PerimPerUnder001=sum(g.per001)/sum(g.per);
    total(i).PerimPerUnder0001=sum(g.per0001)/sum(g.per);
    
    var=abun_rshp.ar(rshp_msk(:,i))/1e6; var=var(:);
    pd(i)=fitdist(var, 'GeneralizedPareto', 'Theta', 0.99*minSize/1e6);
    total(i).a=pd(i).sigma; % size param
    total(i).c=pd(i).k; % shape param
    total(i).k=pd(i).theta;
    par_cdf{i}=gpcdf(ev_ar,pd(i).k,pd(i).sigma,pd(i).theta);
    par_pdf{i}=gppdf(ev_ar,pd(i).k,pd(i).sigma,pd(i).theta);
%     par_cdf_per{i}=gpcdf(ev_per,pd(i).k,pd(i).sigma,pd(i).theta);
%     par_pdf_per{i}=gppdf(ev_per,pd(i).k,pd(i).sigma,pd(i).theta);    
    % plot histogram of area stats
    %       figure
    %     [N,edges] = histcounts(X,edges)
    
        % area
    if plotArea
    figure%(1)
        histogram(abun_rshp.ar(rshp_msk(:,i))/1e6, ev_ar, 'FaceColor','auto', 'Normalization', 'count'); xlabel('Area ($km^2$)'); ylabel('Count');
        title({'Area distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
        set(gca, 'YScale', 'log', 'XScale', 'log')
        annotation(gcf,'textbox',...
            [0.72 0.70 0.25 0.16],...
            'String',['n = ', num2str(sum(rshp_msk(:,i)))],...
            'LineStyle','none',...
            'FontSize',19,...
            'FitBoxToText','off');
        xlim([0 10])
        if plotFit
        hold on
            plot(ev_ar, par_pdf{i})
            legend({labels{i}, 'Pareto PDF fit'}, 'location', 'best')
        hold off
        end
    end
    
    if plotCumArea % only make these plots for total extent  
            % cumulative area
                    figure(2); hold on
            [N,edg]=histcounts(abun_rshp.ar(rshp_msk(:,i))/1e6, ev_ar, ...
                'Normalization', 'cumcount');
            plot(edg(1:end-1), max(N)-N); xlabel('Area ($km^2$)'); ylabel('Count of lakes greater than given area');
            
            
%         figure(2); hold on
            

            
                % plot pareto fit
            if plotFit
                hold off
                [CdfY,CdfX] = ecdf(abun_rshp.ar(rshp_msk(:,i))/1e6,'Function','survivor'); 
                plot(CdfX, CdfY)
                hold on
                YPlot = cdf(pd(i),ev_ar);
                YPlot = 1 - YPlot;
                hLine = plot(ev_ar,YPlot)
                legend({'AirSWOT extent', 'Pareto CDF fit'}, 'location', 'best')
            end
            hold off
            xlabel('Area ($km^2$)'); ylabel('Number of lakes of greater area');
            title({['Cumulative Area distribution (lakes ', num2str(minSize),' - ',num2str(maxSize),' $km^2$)']})%, ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            if i==i_end % last time
                if plotRegions
                    legend(labels, 'location', 'best', 'FontSize', 15)
                else
                end
            end
    end
    
    if plotPerim
            % perim
        figure%(3)
            h=histogram(abun_rshp.per(rshp_msk(:,i))/1e3, ev_ar, 'FaceColor','auto'); xlabel('Perimeter (km)'); ylabel('Count');
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
            [N,edg]=histcounts(abun_rshp.per(rshp_msk(:,i))/1e3, ev_per, ...
                'Normalization', 'cumcount');
            plot(edg(1:end-1), max(N)-N); xlabel('Perimeter km)'); ylabel('Count of lakes greater than given perimeter');
            
            
%         figure(2); hold on
            

            
                % plot pareto fit
            if plotFit
                hold off
                [CdfY,CdfX] = ecdf(abun_rshp.per(rshp_msk(:,i))/1e3,'Function','survivor'); 
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
                'String',['n = ', num2str(sum(rshp_msk(:,i)))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
    end
    
    if plotAreaDevel
                % area vs SDI
        figure%(5)
            plot(abun_rshp.ar(rshp_msk(:,i))/1e6, abun_rshp.SDF(rshp_msk(:,i)), '.')
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

if prelimPlots % additional summary plots
    figure
    histogram([abun.freq_min40]); title('Distribution of lake densities for ABoVE tiles')
    xlabel('Lakes per 100 $km^2$')
    ylabel('Count')
    
    figure
    histogram([abun.lim]); title('Distribution of water fractions')
    xlabel('Limnicity per ABoVE tile (\%)')
    ylabel('Count')
    
    figure
    plot(abun_rshp.SDF, 1./abun_rshp.LehnerDevel, '.')
    xlabel('SDF'); ylabel('Lehner SDF')
    
    figure
    bar([total.lim]*100)
    set(gca, 'XTickLabel', {total.region}, 'XTickLabelRotation', 45)
    title('Water fraction by region')
end

if plotMacDonald %( 3 is Yukon flats)
    i=3;
    % here
    ev_mac=logspace(1.5, 7, 40);
    histogram(abun_rshp.ar(rshp_msk(:,i)), ev_mac, 'FaceColor','auto', 'Normalization', 'count'); xlabel('Area ($m^2$)'); ylabel('Count');
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
    plotLogspace=1;
        % params
    if plotLogspace
%         edges1_area=[4e-5 1e-4 1e-3 0.01 0.1 1];
%         edges1_area=[ min(abun_rshp.ar/1e6) 4e-3 0.04 0.4 max(abun_rshp.ar/1e6)+1];
        edges1_area=logspace(log10(min(abun_rshp.ar/1e6)), 1, 50);
%         edges1_perim=[min(abun_rshp.per/1e3) 0.19 1.9 19  max(abun_rshp.per/1e3)+1];
        edges1_perim=logspace(log10(min(abun_rshp.per/1e3)), 20, 50);        
    else
        edges1_area=linspace((min(abun_rshp.ar/1e6)), 1, 50);
        edges1_perim=linspace((min(abun_rshp.per/1e3)), 20, 50);
    end
    
    
        % just area
    
    figure
    subplot(211)
    [counts, ~, binsArea] = histcounts(abun_rshp.ar/1e6,edges1_area);
    histogram(abun_rshp.ar/1e6,edges1_area);
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    end
    xlabel('Area ($km^2$)'); ylabel('count')
    title('Area histogram')
    
    subplot(212)
    aggAreabyArea=zeros(size(counts)); % init
    for i=1:max(binsArea)
        aggAreabyArea(i)=sum(abun_rshp.ar(binsArea==i))/1e6;
        
    end
    if sum(binsArea==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsArea==0))
    end
    histogram('BinEdges',edges1_area, 'BinCounts', aggAreabyArea)
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    end
    xlabel('Area ($km^2$)'); ylabel('Sum of binned areas ($km^2$)')
    title('Area histogram binned by area') 
    
        % just perim
    figure
    subplot(211)
    histogram(abun_rshp.per/1e3,edges1_perim);
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    end
    xlabel('Perimeter ($km$)'); ylabel('count')
    title('Perimeter histogram')
        
    subplot(212)
    [counts, ~, binsPerim] = histcounts(abun_rshp.per/1e3,edges1_perim);
    aggPerimbyPerim=zeros(size(counts)); % init
    for i=1:max(binsPerim)
        aggPerimbyPerim(i)=sum(abun_rshp.per(binsPerim==i))/1000;
        
    end
    if sum(binsPerim==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsPerim==0))
    end
    histogram('BinEdges',edges1_perim, 'BinCounts', aggPerimbyPerim)
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    end
    xlabel('Perim ($km^2$)'); ylabel('Sum of binned perimeters (km)')
    title('Perimeter histogram binned by perimeter')
    
        % bin perim by area
    figure    
    aggPerimbyArea=zeros(size(counts)); % init
    for i=1:max(binsArea)
        aggPerimbyArea(i)=sum(abun_rshp.per(binsArea==i))/1000;
        
    end
    if sum(binsArea==0)>0
            warning('Not all data was binned.')
            fprintf('\tFigure: %d, i= %d\n', get(gcf,'Number'), i)
            disp(sum(binsPerim==0))
        end
    histogram('BinEdges',edges1_area, 'BinCounts', aggPerimbyArea)
    if plotLogspace
        set(gca, 'YScale', 'lin', 'XScale', 'log')
    end
    xlabel('Area ($km^2$)'); ylabel('Sum of binned perimeters (km)')
    title('Perimeter histogram binned by area')
end

%% save figs
if saveFigs
for i=1:get(gcf, 'Number')
    saveas(i, [figs_out, 'Geom_', num2str(i), '.png'])
end

    % save text file pointing to current directory (for this script)
fid=fopen([figs_out, 'Source.txt'], 'w+');
fprintf(fid, '%s\n', pwd);
fclose(fid)
end
% strrep(total(i).region, ' ',''),


%% output shapefile

if saveShp
    disp('Saving shape...')
    S=mappoint(abun_rshp.long,abun_rshp.lat, abun_rshp);
    shapewrite(S, shp_out);
end


%% computations


%% output stats table

% tbl=struct2table(total);
% writetable(tbl, tbl_out);