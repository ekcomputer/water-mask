function make_PL_plot(list,alpha,xmin,num_per_decade,pvals,varargin)
% uses 0.1 as boundary bw possible and impossible PL
if ~isempty(varargin)
    
    do_ebars = 1;
    ebars = varargin{1};
    
end


horvat_colors;
lw=3; % line width
minA = xmin(1);
maxA = max(list);
Alist = list(list > minA);
maxlog = ceil(log10(max(list)));
minlog = floor(log10(min(list)));
width = maxlog - minlog;
    % EK added
stem_ebars=0; % 1 = stem plot; 2=bar plot

Abins = logspace(log10(minA),log10(minA)+width,num_per_decade*width);
% Center of those bins
Acenters = Abins(1:end-1) + diff(Abins)/2;

% Take the histogram of those floes. If there are none, don't display as
% nan but make small on the plot.
to_plot = histcounts(Alist,Abins);
to_plot(to_plot == 0) = .01;

% Find the location from which to draw the tail lines - closest to
[~,b] = min(abs(Acenters - xmin(2)));

% Make the pure-power-law distributions

% For all floes via the V-C method
distro_all = diff(Abins).*Acenters.^(-alpha(1));
distro_all = to_plot(1) .* distro_all./distro_all(1);

% For the tail for the V-C method
distro_mle_tail = diff(Abins(b:end)).*Acenters(b:end).^(-alpha(2));
distro_mle_tail = to_plot(b).*distro_mle_tail./distro_mle_tail(1);

if do_ebars
    
    distro_e(:,1) = diff(Abins(b:end)).*Acenters(b:end).^(-alpha(2) + ebars(1));
    distro_e(:,1) = to_plot(b).*distro_e(:,1)./distro_e(1,1);
    distro_e(:,2) = diff(Abins(b:end)).*Acenters(b:end).^(-alpha(2) - ebars(1));
    distro_e(:,2) = to_plot(b).*distro_e(:,2)./distro_e(1,2);
    
end

% Difference at the tail

diff_all = distro_all - to_plot(1:end);
diff_tail = distro_mle_tail - to_plot(b:end);

misfit_tail = round(100*sum(abs(diff_all(2:end)) / sum(distro_all(2:end))));
fprintf('The full misfit is %d percent \n ',misfit_tail);

misfit_tail = round(100*sum(abs(diff_tail(2:end)) / sum(distro_mle_tail(2:end))));
fprintf('The tail misfit is %d percent \n ',misfit_tail);

% Plot the beginning of the tail.

% The actual histogram
p_data = loglog(Acenters,to_plot,'linewidth',lw,'color','k');
hold on

    % add uncertainty in red line (X0)
if do_ebars & pvals(2)>=0.1
    xvals = logspace(log10(Acenters(b) - ebars(2)),log10(Acenters(b) + ebars(2)),10);  
    area(xvals,10.^ceil(log10(max(to_plot))) + 0*xvals,'facecolor',.9*[1 1 1],'facealpha',1,'linestyle','none')
end

p_mle = loglog(Acenters,distro_all,'linewidth',lw,'color',clabs(2,:));
p_tail = loglog(Acenters(b:end),distro_mle_tail,'linewidth',lw,'color',clabs(1,:));

if do_ebars
    loglog(Acenters(b:end),distro_e,'--','linewidth',0.5*lw,'color',0.5*clabs(1,:)); 
end

if stem_ebars==1
    stem(Acenters,abs(diff_all), 'marker','+','markersize', 8,'markerfacecolor',clabs(2,:), 'LineWidth', 4);
    stem(Acenters(b:end),abs(diff_tail), 'marker','+','markersize', 8,'markerfacecolor',clabs(1,:), 'LineWidth', 4);
elseif stem_ebars==2
    bar(Acenters,abs(diff_all),'facecolor',clabs(2,:));
    bar(Acenters(b:end),abs(diff_tail), 'facecolor',clabs(1,:));
else
    area(Acenters,abs(diff_all),'facealpha',0.25,'facecolor',clabs(2,:),'linestyle','none');
    area(Acenters(b:end),abs(diff_tail),'facealpha',0.5,'facecolor',clabs(1,:),'linestyle','none');
end
%     scatter(Acenters,0*Acenters + 1,'+k')

    % add vertical red line (X0)
if pvals(2)>=0.1
    plot(0*Acenters + Acenters(b),logspace(0,7,length(Acenters)),'color','r');
end
h = get(gca,'children');
h = [h(end); h(1:end-1)];
set(gca,'children',h);

ylim([1 10.^ceil(log10(max(to_plot)))])
xlim([minA maxA]);
grid on
box on
set(gca,'xscale','log','yscale','log');

legend([p_data p_mle p_tail],{'Data',sprintf('p=%1.2f',pvals(1)),sprintf('p=%1.2f',pvals(2))})