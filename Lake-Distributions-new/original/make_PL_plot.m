function make_PL_plot(list,alpha,xmin,num_per_decade,pvals,varargin)

if ~isempty(varargin)
    
    do_ebars = 1;
    ebars = varargin{1};
    
end


horvat_colors;

minA = xmin(1);
maxA = max(list);
Alist = list(list > minA);
maxlog = ceil(log10(max(list)));
minlog = floor(log10(min(list)));
width = maxlog - minlog;
    % EK added
stem_ebars=true;

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
p_data = loglog(Acenters,to_plot,'linewidth',1,'color','k');
hold on

if do_ebars
    xvals = logspace(log10(Acenters(b) - ebars(2)),log10(Acenters(b) + ebars(2)),10);  
    area(xvals,10.^ceil(log10(max(to_plot))) + 0*xvals,'facecolor',.9*[1 1 1],'facealpha',1,'linestyle','none')
end

p_mle = loglog(Acenters,distro_all,'linewidth',1,'color',clabs(2,:));
p_tail = loglog(Acenters(b:end),distro_mle_tail,'linewidth',1,'color',clabs(1,:));

if do_ebars
    loglog(Acenters(b:end),distro_e,'--','linewidth',0.5,'color',0.5*clabs(1,:)); 
end

area(Acenters,abs(diff_all),'facealpha',0.25,'facecolor',clabs(2,:),'linestyle','none');
area(Acenters(b:end),abs(diff_tail),'facealpha',0.5,'facecolor',clabs(1,:),'linestyle','none');
scatter(Acenters,0*Acenters + 1,'+k')

plot(0*Acenters + Acenters(b),logspace(0,7,length(Acenters)),'color','r');

h = get(gca,'children');
h = [h(end); h(1:end-1)];
set(gca,'children',h);

ylim([1 10.^ceil(log10(max(to_plot)))])
xlim([minA maxA]);
grid on
box on
set(gca,'xscale','log','yscale','log');

legend([p_data p_mle p_tail],{'Data',sprintf('p=%d',round(pvals(1)*100)),sprintf('p=%d',round(pvals(2)*100))})