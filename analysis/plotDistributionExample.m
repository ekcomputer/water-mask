% plots distributions for pareto, lognormal, etc

    % set x - doesn't matter
% ev=logspace(-4, 0, 200);
% ev=logspace(log10(0.016), 0, 200);
ev=linspace(0, 1, 200);
% p=cdf('Generalized Pareto', ev);

    % gen fxns and plot
clf; hold on
p=pdf('Generalized Pareto', ev, 1,1,0.02);
l=pdf('Lognormal', ev, 0.0001, 1.5);
poi=pdf('Poisson', ev, 1);
e=pdf('Exponential', ev, 1);
g=pdf('Gamma', ev, 1, 1);
plot(ev,p, 'LineWidth', 5)
plot(ev, l, 'LineWidth', 5)
% plot(ev, poi)
plot(ev, e, 'LineWidth', 5)
% plot(ev, g)
hold off


    % view params
set(gca, 'YScale', 'lin', 'XScale', 'lin')
% text(max(ev),p(end),'label', 'FontSize', 18)