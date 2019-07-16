% helper to fix problem of file getting too big , mid-run
% just used once during a run with a bug

for i=1:length(abun)
    abun(i).stats=rmfield(abun(i).stats, 'PixelList');
end
%%
for i=1:length(abun)
    try
        e(i)=~isempty(abun(i).stats(1).SDF);
    catch
        e(i)=0;
        [abun(i).stats.SDF, abun(i).stats.LehnerDevel] = [];
    end
end

%%
for i=1:length(abun)
    fprintf('%d   %d\n', i,length(abun(i).stats) )
    if i==39 || i==53 || i==54
        abun(i).stats=ex;
    end
end