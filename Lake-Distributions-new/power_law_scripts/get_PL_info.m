function [alpha,xmin,LL,pval,gof,ebar] = get_PL_info(Areas)

alpha = [0 0];
pval = [0 0];
gof = [0 0];

% repnum = 250;  % set to 1000, will increase time by 64x
% samp_num = 250; 
% semi_param_number = 250; 
global env
repnum = env.plinfo.repnum;  % set to 1000, will increase time by 64x
samp_num = env.plinfo.samp_num; 
semi_param_number = env.plinfo.semi_param_number; 

% Find the peak of the distribution
xmin(1) = get_top_of_distribution(Areas);
Areas = Areas(Areas >= xmin(1)); 
% Operate on all lake areas in that bin
% Q1 - Look at all the lakes areas. If we include all, how good is the fit?
% Do we get a full spectrum of power laws
[alpha(1),~,~] = plfit(Areas,'xmin',xmin(1));
fprintf('Distribution peak is at %.02f: decay coefficient above that is %.02f. \n',xmin(1),alpha(1));

% Q2 - ok what about just the tail? We let the pl code handle
% determining the tail.
[alpha(2),xmin(2),LL] = plfit(Areas);
[ebar(1),ebar(2),ebar(3)] = plvar(Areas,'reps',repnum,'sample',samp_num); % remove these calls to repnum and samp_num for acc, not speed

above = 100*sum(Areas > xmin(2)) / numel(Areas); 

fprintf('MLE Tail begins at %.02f: decay coefficient above that is %.02f. Accounts for %.02f percent of DS \n',xmin(2),alpha(2),above);

% Q3 - how good is this fit for all lakes
% [pval(1),gof(1)]=plpva(Areas,xmin(1),'reps',repnum,'sample',semi_param_number,'silent');
[pval(1),gof(1)]=plpva(Areas,xmin(1),'reps',repnum,'xmin',xmin(1),'sample',semi_param_number); % remove semi_param

fprintf('For all, p val is %.02f. ',pval(1));

% Q4 what about for the tail?

trunc_A = Areas(Areas > xmin(2));

% [pval(2),gof(2)]= plpva(Areas,xmin(2),'reps',repnum,'sample',semi_param_number,'silent');
[pval(2),gof(2)]= plpva(trunc_A,xmin(2),'reps',repnum,'sample',semi_param_number,'xmin',xmin(2)); % remove semi_param

fprintf('For the tail, p val is %.02f \n',pval(2));
