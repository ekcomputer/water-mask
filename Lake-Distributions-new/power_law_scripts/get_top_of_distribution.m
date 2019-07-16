function b_min_all = get_top_of_distribution(numlist)

%%
mult = 1/min(numlist); 

numlist = mult*numlist; 

bins = 0:1:50;
numlist = round(numlist);   
num = histc(numlist,bins);
[~,ind] = max(num);
b_min_all = bins(ind)/mult;


end