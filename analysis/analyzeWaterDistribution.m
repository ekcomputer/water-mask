% function to load and analyze water distribution in aggregate and by
% region.  Also creates shapefile

%% params
clear
saveFigs=0;
%% directories
if ~isunix
    struct_in='J:\output\analysis\distrib.mat';
else
    struct_in='/Volumes/Galadriel/output/analysis/distrib.mat';
    figs_out='/Volumes/Galadriel/output/pic/geomFigsBulk/';
end

%% load and reshape input data
load(struct_in);
abun_rshp.ar=[]; abun_rshp.per=[]; abun_rshp.lat=[]; abun_rshp.file_idx=[];
for i=1:length(abun)
    ar_temp=[abun(i).stats.Area];
    per_temp=[abun(i).stats.Perimeter];
    lat_temp=[abun(i).stats.lat];
    file_idx_temp=repmat([abun(i).file_idx], [1,length(abun(i).stats)]);
    abun_rshp.ar=[abun_rshp.ar,ar_temp];
    abun_rshp.per=[abun_rshp.per,per_temp];
    abun_rshp.lat=[abun_rshp.lat,lat_temp];
    abun_rshp.file_idx=[abun_rshp.file_idx,file_idx_temp]; % record index for masking
end

%% regions

% all
regions{1}=1:330;
labels{1}='AirSWOT extent';

% north slope
regions{2}=269:272;
labels{2}='North Slope';

% yukon flats
regions{3}=1+[65	66	67	68	69	70	71	72	73	74	90	91	263	264	265	266	267	268	273	274	275	276	277	278	279	280];
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

% regions{8}=1+
% regions{9}=1+
% regions{10}=1+
% regions{11}=1+
% regions{3}=1+
% regions{3}=1+
% regions{3}=1+
% regions{3}=1+

%% plot
close all
makePlot=[1 1 1 1 1]; % which plots to draw
% ev=sort(10-logspace(-4, 1, 200)); % edge vector to use for binning
ev=logspace(-4, 1, 200);
for i=1:length(regions)
    regions{i}=intersect(regions{i}, [abun.file_idx]); % use only idx with observations included in complete map
    msk=false(size(abun_rshp.file_idx));
    for j=1:length(regions{i})
        msk=msk | abun_rshp.file_idx==regions{i}(j); % true values are regions (files) of interet
    end
    % plot histogram of area stats
    %       figure
    %     [N,edges] = histcounts(X,edges)
    
        % area
    if makePlot(1)
    figure%(1)
        histogram(abun_rshp.ar(msk)/1e6, ev, 'FaceColor','auto', 'Normalization', 'probability'); xlabel('Area ($km^2$)'); ylabel('Count');
        title({'Area distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
        set(gca, 'YScale', 'log', 'XScale', 'log')
        annotation(gcf,'textbox',...
            [0.72 0.70 0.25 0.16],...
            'String',['n = ', num2str(sum(msk))],...
            'LineStyle','none',...
            'FontSize',19,...
            'FitBoxToText','off');
        xlim([0 10])
    end
    if makePlot(2) % only make these plots for total extent  
            % cumulative area
        figure%(2)
            [N,edg]=histcounts(abun_rshp.ar(msk)/1e6, ev, ...
                'Normalization', 'cumcount');
            plot(edg(1:end-1), max(N)-N); xlabel('Area ($km^2$)'); ylabel('Count of lakes greater than given area');
            xlabel('Area ($km^2$)'); ylabel('Count'); title({'Cumulative Area distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(sum(msk))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
    end
    if makePlot(3)
            % perim
        figure%(3)
            h=histogram(abun_rshp.per(msk)/1e3, ev, 'FaceColor','auto'); xlabel('Perimeter (km)'); ylabel('Count');
            title({'Perimeter distribution', ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(sum(msk))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
            xlim([0.01 10])
    end
    if makePlot(4)
                % perim/area
        figure%(4)
            h=histogram(abun_rshp.per(msk)/1e3./(abun_rshp.ar(msk)/1e6), 'FaceColor','auto'); xlabel('1/Length (1/km)'); ylabel('Count');
            title({'Perimeter:area distribution' , ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'linear', 'XScale', 'linear')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(sum(msk))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
    end
    
    if makePlot(5)
                % area vs perim/area
        figure%(5)
            plot(abun_rshp.ar(msk)/1e6, (abun_rshp.per(msk)/1e6)./(abun_rshp.ar(msk)/1e6), '.')
            ylabel('1/Length (1/km)'); xlabel('Area ($km^2$)');
            title({'Perimeter:area distribution' , ['region: ', labels{i}]}, 'Interpreter', 'none')
            set(gca, 'YScale', 'log', 'XScale', 'log')
            annotation(gcf,'textbox',...
                [0.72 0.70 0.25 0.16],...
                'String',['n = ', num2str(sum(msk))],...
                'LineStyle','none',...
                'FontSize',19,...
                'FitBoxToText','off');
    end
end

%% save figs
if saveFigs
for i=1:get(gcf, 'Number')
    saveas(i, [figs_out, 'GeomFig_', num2str(i), '.png'])
end

    % save text file pointing to current directory (for this script)
fid=fopen([figs_out, 'Source.txt'], 'w+');
fprintf(fid, '%s\n', pwd);
end
%% count total land and water

total.land=sum([abun.land])/1e6;
total.water=sum([abun.water])/1e6;
total.area=total.land + total.water;
total.lim=total.water/(total.water+total.land);
total.lim2=mean([abun.lim].*([abun.land]+[abun.water])/sum([abun.land]+[abun.water])) % double check...




