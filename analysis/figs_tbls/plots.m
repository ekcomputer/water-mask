% script to make plots and tables for paper

%% load data
clear; close all
global plt env

saved_dir=env.saved_dir;
tbl_out=env.tbl_out;
geom_in=env.geom_in;
labels_in=env.labels_in;
labels_exp_in=env.labels_exp_in;

load(labels_in);
load(labels_exp_in);
plt.c=[0 0 0.8];

plot_morph=1;
plot_gcp=1;
plot_pl_all=0;
plot_pl_region=1;
alignPlots=1;
plotHist=1;
calcs=0;
makeTable=0;
lw=1.5; % line width
rot=0; % label rotation
regions_Q=env.regions_Q; % for power law plots and table
cat_Q=env.cat_Q;%[1:15]; % labels(Q) to plot

    % load .mat vars
% files=cellstr(ls([saved_dir, '*.mat']));
% for i=1:length(files)
%     load([saved_dir,files{i}])
% end
load(env.geom_in);
load(env.analyzeWaterDistribution); load(env.fit_data_reg);
% load(env.fit_data); 
load(env.plotGCPAcc); load(env.labels_exp_in);
% labels{24}={'Pothole'; 'Lakes'};
% labels{22}='Pothole \newline Lakes';
% labels{23}='Shield \newline Lakes';
% labels{24}='Wetland- \newline Lakes';
% labels{25}='Thermokarst \newline Lakes';
% labels{25}='Valley \newline Lakes';
total%% morphometry plots
if plot_morph

    figure
    b(1)=bar([total(cat_Q).perUnder001]*100, 'FaceColor', 'flat');
    set(gca, 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', rot)
    title('Percent of water bodies under 0.001 $km^2$')
    ylabel('Percent')
%     b.CData(2,:) = [0 0.8 0.8];
    
    figure
    b(2)=bar([total(cat_Q).ArPerUnder001]*100, 'FaceColor', 'flat');
    set(gca, 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', rot)
    title('Percent of areas under 0.001 $km^2$')
    ylabel('Percent')
    
    figure
    b(3)=bar([total(cat_Q).PerimPerUnder001]*100, 'FaceColor', 'flat');
    set(gca, 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', rot)
    title('Perimeters from water bodies under 0.001 $km^2$')
    ylabel('Percent')
  
    figure
    b(4)=bar([total(cat_Q).MedArea], 'FaceColor', 'flat');
    set(gca, 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', rot)
    title('Median water body area')
    ylabel('$km^2$')

    
    figure
    b(5)=bar([total(cat_Q).MedPerim], 'FaceColor', 'flat');
    set(gca, 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', rot)
    title('Median perimeter')
    ylabel('($km$)')
    
    figure;
    b(6)=scatter(g(1).ar, SDF_all);
    xlabel('Area ($km^2$)'); ylabel('SDF')
    set(gca, 'YScale', 'log', 'XScale', 'log')
    
    figure
    b(7)=bar([total(cat_Q).MedSDF], 'FaceColor', 'flat');
    ylabel('SDF')
    set(gca, 'XTickLabel', labels(cat_Q), 'XTick', 1:length(cat_Q), 'XTickLabelRotation', rot)
    title('Median SDF')
    
    figure
    x=[total(regions_Q).MedArea];
    y=[total(regions_Q).MedSDF];
    b(8)=scatter(x, y, 'LineWidth', 2);
    xlabel('Median area ($km^2$)'); ylabel('SDF')
        % add text labels
    dx=0.00006; dy=0.006; % offsets for labels
    text(x+dx, y+dy,...
        labels(regions_Q), 'Interpreter', 'latex', 'FontSize', 19);
        % curve fitting: use power law
    cfit=fit(x(:),y(:),'power1');
    fit_func=@(x) cfit.a.*x.^(cfit.b);
    xs=linspace(min(x), max(x), 100);
    hold on; plot(xs,(fit_func(xs)),':', 'Color', [0.8 0.8 0.8])
    box on
    
    figure
    vals=[geom(cat_Q).fraction_water]*100;
    b(9)=bar(vals); b(9).FaceColor=plt.c;
    title('Open water fraction')
    labs=labels(cat_Q); labs{3}={'Pothole', 'Lakes'};
    for i=1:numel(labs)
    text(i,vals(i),num2str(vals(i),'%0.2f'),...
               'HorizontalAlignment','center',...
               'VerticalAlignment','bottom', 'fontsize', 19)
    end
    ylabel('Percent')
    set(gca, 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', rot)
    set(gcf, 'Position', [-1344 -379 1113 784])
    
end    
    %% GCP
if plot_gcp
    err=gcp.data(:,8);
    figure;
    b(10)=histogram(err); b(10).FaceColor=plt.c; b(10).LineWidth=lw;
    b(10).FaceAlpha=1
    xlabel('Geolocation error (m)')
    ylabel('Count')
    set(gca, 'LineWidth', 1)
end
%% plot hist

if plotHist
        % load any figs
    files=cellstr(ls([saved_dir, '*.fig']));
    for i=1:length(files)
        H=open([saved_dir,files{i}])
    end
        % touch up
        xlabel('NDWI (DN)'); ylabel('Number of pixels')
        title('')
        set(gcf, 'WindowStyle', 'normal')
end

%% touch up my figs
if plot_morph || plotHist
    for d=1:get(gcf, 'Number') %h =  findobj('type','figure'); n = length(h);
       set(gca, 'TitleFontSizeMultiplier', 1.5)
        figure(d)
       box on; 
       if d~=8
           grid on
       end
       set(gca, 'LineWidth', 1, 'FontSize', 20) % from 25
       set(gcf, 'Position', [-1344        -341         767         746]) % fix label
       child=get(gca, 'Children'); 
       if ~strcmp(child(end).Type, 'bar')
       else
          xlim([0.5,length(cat_Q)+0.5])
          set(gca, 'XTick', [1:length(cat_Q)]);
       end
                set(gca, 'LineWidth', lw, 'YMinorTick', 'on',...
                    'GridAlpha', 0.25, 'TitleFontSizeMultiplier', 1, 'TitleFontWeight', 'bold',...
                    'FontName', 'Ariel');
       if d~=6 && d~=8 && d~=9
           try b(d).CData=repmat(plt.c, [size(b(d).CData, 1),1, 1]); end
           if cat_Q(1)==1 % if showing all
               b(d).CData(1,:)=[0.71 1 1];
           end
           try b(d).LineWidth=lw+0.5; end
       else
           try
                b(d).CData=plt.c;
           catch
               b(d).FaceColor=plt.c;
           end
       end
    end
end

%% power law plots
numperdecade = 3; 
if plot_pl_all
    addpath('D:\Dropbox\Matlab\Lake-Distributions-new\power_law_scripts\');
    addpath('D:\Dropbox\Matlab\Lake-Distributions-new\');

    figure
%     make_PL_plot(Fused_area,Fused_alpha,Fused_xmin,numperdecade,Fused_pval,Fused_ebar)
    make_PL_plot_FINESST(Fused_regional{i},alpha_regional(1,:),xmin_regional(1,:),numperdecade,pval_regional(1,:),ebar_regional(1,:))
    xlabel('Area ($km^2$)'); ylabel('Count'); 
%     title(sprintf('\\textbf{%s}\nn = %d, $\\alpha$ = %0.2f','All fused water bodies',round(Fused_ebar(3)), Fused_alpha(2)));
    % title(sprintf('Fused Data: n = %d',round(Fused_ebar(3))));
    set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100 1000], 'FontSize', 21,...
        'LineWidth', 1.5)
    pos = [12 8]; pos=[10 8];
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');
    set(gcf,'windowstyle','normal','position',[0 0 pos],'paperposition',[0 0 pos],'papersize',pos,'units','inches','paperunits','inches');

        % multiple

    % set(0,'DefaultLineMarkerSize',6);
    % set(groot,'defaultAxesFontSize',11); %EK I changed this to 14
    set(groot,'defaultLineLineWidth',1);
    grid minor % turn off minor ticks
    axis ([1e-4 1e3 1 1e4])
end
   pl_regions={regions_Q, [cat_Q, 1], 1};
if plot_pl_region
        % subplot params (gap_h, gap_w)
    p(1).gap=[0.057, 0.021]; p(1).marg_h=0.06; p(1).marg_w=0.03; p(1).fontsize=14;
    p(2).gap=[0.07, 0.029]; p(2).marg_h=0.07; p(2).marg_w=0.04; p(2).fontsize=17;
    for j=1:2%length(pl_regions)
        figure
        c=1;
        for i = pl_regions{j}
            if c==1
            if length(pl_regions{j}) <10 % just categories
%                 subplot(2,3,c)
                ha=tight_subplot(2,3,p(j).gap,p(j).marg_h,p(j).marg_w)              
            else % all regions
%                 subplot(3,5,c)
                ha=tight_subplot(3,5,p(j).gap,p(j).marg_h,p(j).marg_w)
                delete(ha(14)); delete(ha(15))
            end
            end
            set(gcf,'CurrentAxes',ha(c))
            make_PL_plot(Fused_regional{i},alpha_regional(i,:),xmin_regional(i,:),numperdecade,pval_regional(i,:),ebar_regional(i,:)) %ebar_regional(i,:)
            %     xlabel('Area ($km^2$)'); ylabel('Count'); 
            if j==1
                title(sprintf('\\fontsize{12}{18}\\textbf{%s}\nn = %d, $\\alpha$ = %0.2f',...
                    labels{(i)},sum(Fused_regional{i}>=xmin_regional(i,2)), alpha_regional(i,2)));
            else
                title(sprintf('\\fontsize{12}{18}\\textbf{%s}\nn = %d, $\\alpha$ = %0.2f',...
                    labels_expl{(i)},sum(Fused_regional{i}>=xmin_regional(i,2)), alpha_regional(i,2)));
            end
            set(gca, 'XTick', [0.0001 0.001 0.01 0.1 1 10 100 1000], 'LineWidth', 1.5,...
                'FontSize', p(j).fontsize, 'TitleFontSizeMultiplier', 1)
%                  'XTickLabel', {'-4','-2','-2','-1','0','1','2','3'})
            grid minor % turn off minor ticks
            axis ([1e-4 1e3 1 1e4])
            if (j==1 && ~ismember(c, [1,6,11])) || (j==2 && ~ismember(c,[1,4]))
                set(gca, 'YTickLabel', '')
            end
            if (j==1 && ~ismember(c, [9 10 11 12 13])) || (j==2 && ~ismember(c,[4 5 6]))
                set(gca, 'XTickLabel', '')
            end
            c=c+1;
        end


        pos = {[1441 875 ], [ 1280 905]}; %,[ 1108 758],[ 1127 778 ]}; % set figure aspect
            % old position: [1178 1018]
        set(gcf,'windowstyle','normal','position',[ -1448.5       -524.5     pos{j}])
    end
end
set(groot,'defaultLineLineWidth',4);
%% move onto axes in alignment
if alignPlots && plot_morph
    subQ={[1 2 4 9]};
    p(3).gap=[0.06, 0.052]; p(3).marg_h=[0.13 0.06]; p(3).marg_w=[0.05 0.02];
    for j=1:length(subQ)
        figure;
        c=1; %counter
        for i=subQ{j}
            if c==1
                sub=tight_subplot(2,2,p(3).gap,p(3).marg_h,p(3).marg_w) 
%                 sub(c)=subplot(2, 2,c);
                set(gcf, 'Position', [-1428        -509        1251         944])
            end
            sub(c).TitleFontSizeMultiplier=1.5; %%% <------------- HERE wrong shape
            copyobj(b(i),sub(c))
            set(sub(c), 'FontSize', 12, 'XLabel', get(b(i).Parent, 'XLabel'), 'YLabel', get(b(i).Parent, 'YLabel'),...
                'Title', get(b(i).Parent, 'Title'));
                % font size for title, indep. of paretn
            set(sub(c), 'TitleFontSizeMultiplier', 1.5);
            set(b(i).Parent, 'XLabel', get(sub(c), 'XLabel'),...
                'Title', get(sub(c), 'Title'));
            set(gcf,'CurrentAxes',sub(c))
            xlim([0.5,length(cat_Q)+0.5])
            box on
            set(sub(c), 'LineWidth', lw, 'XTick', [1:length(cat_Q)], 'YMinorTick', 'on',...
                'GridAlpha', 0.25, 'TitleFontSizeMultiplier', 1, 'TitleFontWeight', 'bold',...
                'FontName', 'Ariel');
            grid on
            set(sub(c), 'XTickLabel', labels(cat_Q), 'XTickLabelRotation', 40,...
                'FontSize', 20);
            set(sub(c), 'TitleFontSizeMultiplier',...
                1.3, 'LabelFontSizeMultiplier', 1.5);
%             if c==length(subQ{j}) % if end
%                 set(gca, 'XTickLabel', labels(Q), 'XTickLabelRotation', 0);
%             else
%                 set(gca, 'XTickLabel', {});
%             end
%             if ismember(c, [2,4])
%                 set(sub(c), 'YTickLabel', '')
%             end
            if ~ismember(c, [3,4])
                set(sub(c), 'XTickLabel', '')
            end
            c=c+1;
        end
        
    end
end
set(gcf,'windowstyle','normal','position',[ -1448.5   -524.5  1198 1031])

%% calculate % covered by PLH, with errors propogated
if calcs
        % close figures that were consolidated
    for k=[1 2 3 4 5 7 9]; close(figure(k)); end
    
    load(env.lake_databases, 'hl_fused')
    for j=1:3 % lower, value, and upper bounds
        for i=[1, env.regions_Q, env.cat_Q]
            if j==1
                a0=xmin_regional(i,2)-ebar_regional(i,2);      
            elseif j==2
                a0=xmin_regional(i,2);
            elseif j==3
                 a0=xmin_regional(i,2)+ebar_regional(i,2);       
            end
            if i==1
                total_area=[hl_fused([hl_fused.Region4]>0).Area];
                d=[hl_fused([hl_fused.Area]>=a0).Area];
            elseif i>20
                total_area=[hl_fused([hl_fused.Category4]==i-(min(env.cat_Q)-1)).Area];
                d=[hl_fused([hl_fused.Area]>=a0 & [hl_fused.Region4]==i-(min(env.cat_Q)-1)).Area];
            else
                total_area=[hl_fused([hl_fused.Region4]==i).Area];
                d=[hl_fused([hl_fused.Area]>=a0 & [hl_fused.Region4]==i).Area];
            end

            total(i).pl_a_fraction(j)=sum(d)/sum(total_area)*100;
            total(i).pl_c_fraction(j)=numel(d)/numel(total_area)*100;
        end
    end
end
%% output stats table
if makeTable
    load(geom_in)
    delete(tbl_out)
    stable=total; % summ table
    for k=1:length(stable) % deals are unness...
        
            % vars from power law fits
        [stable(k).p_all]=deal(pval_regional(k,1));    
        [stable(k).p]=deal(pval_regional(k,2)); 

        [stable(k).alpha_all]=deal(alpha_regional(k,1)); 
        [stable(k).alpha]=deal(alpha_regional(k,2)); 

        [stable(k).a_min]=deal(xmin_regional(k,2));   

        [stable(k).a_min_error]=deal(ebar_regional(k,2)); 
        [stable(k).alpha_error]=deal(ebar_regional(k,1)); 
        
        stable(k).a_min=stable(k).a_min*1e6; % convert to m2

            % vars from computeWaterFraction.m
        stable(k).area_region=geom(k).area;
        stable(k).area_water=geom(k).area_water;
        stable(k).area_lakes=geom(k).area_lakes;
        stable(k).count_water_chk=geom(k).count_water_chk;
        stable(k).count_lakes_chk=geom(k).count_lakes_chk;
        stable(k).fraction_water=geom(k).fraction_water*100;
        stable(k).fraction_lakes=geom(k).fraction_lakes*100;
        
            % vars from analyzeWaterDistrib- dec to %
        stable(k).perUnder001=stable(k).perUnder001*100;
        stable(k).ArPerUnder001=stable(k).ArPerUnder001*100;
        stable(k).PerimPerUnder001=stable(k).PerimPerUnder001*100;
        stable(k).lim=stable(k).lim;
        stable(k).MedArea=stable(k).MedArea*1e6;
        
        % lim
%         stable(k). % percents to zeros
    end
    vars_Q=[1,regions_Q, cat_Q];
    tbl=struct2table(stable(vars_Q));
    writetable(tbl, tbl_out,'Sheet', 'AllVars');
        % write second table  
    fnames=fieldnames(stable);
    stable_select=rmfield(stable,setdiff(fnames,(fnames([1 2 5 9 12 15 19 27 28 29 31 34 32 33, 35, 40, 41]))));
    tbl=struct2table(stable_select(vars_Q));
    writetable(tbl, tbl_out,'Sheet', 'SelectVars');
        % write third table
    overview=array2table({labels_expl{vars_Q};labels{vars_Q}}');
    overview.area=[geom(vars_Q).area]';
    writetable(overview, tbl_out,'Sheet', 'Overview');
end
% saveallfigs('Figs_regions4', 'mat', 'eps')