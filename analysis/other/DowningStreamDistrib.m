% script to use data from Downing et al 2012 to plot
% they also include length and width

clear; 

n=[28550000 6000000 1260000 264000 55500 11700 2450 515 110 23 5 1]; % count
a=[36500 39200 39600 42500 72800 88100 76400 74300 82600 64900 25400 19800] 

% plot(n,a, 'xk-'); xlabel('Number'); ylabel('Total Area'); title('Downing et al. 2012')
% figure
% plot(a,n, 'xk-'); ylabel('Number'); xlabel('Total Area'); title('Downing et al. 2012')
% 
% bar(a, log(n))
% bar(a./n, n)
% [avg_a_sort, seq]=sort(a./n);
[total_a_sort, seq]=sort(a, 'ascend');
[n_sort, n_seq]=sort(n);

figure(1)
semilogy(a./n, n); ylabel('Number'); xlabel('Average area $(km^2)$'); title('Downing et al. (2012) Rivers')
figure(2)
semilogy(total_a_sort, n(seq)); ylabel('Number'); xlabel('Total area $(km^2)$'); title('Downing et al. (2012) Rivers')
figure(3)
loglog(n_sort, a(n_seq)); xlabel('Number'); ylabel('Total area $(km^2)$'); title('Downing et al. (2012) Rivers')
figure(4)
loglog(a(n_seq)./n(n_seq), n_sort); ylabel('Number'); xlabel('Avg area $(km^2)$'); title('Downing et al. (2012) Rivers')

% 