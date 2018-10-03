function L=fillLabelGaps(L)
% takes labeld matrix L (double) and if there are gaps in label indexes,
% shifts everything else down to fill gaps.
% warning: memory intensive!  Works as long as range(L(:)) aprox equal to
% twice length(unique(L))
disp('   fillLabelGaps')
if range(L(:)) > 3*length(unique(L))
    warning('Range of L is more than 3x the length(unique(L))')
elseif range(L(:)) > 10*length(unique(L))
    error('Range of L is more than 10x the length(unique(L)).  This will consume too much memory.')
end
pixelIndexList = label2idx(L); sz=size(L); clear L
CC=bwconncomp(zeros(sz)); 
c=1; %counter
disp('       loop')
for i=1:length(pixelIndexList)
    if isempty(pixelIndexList{i})
    else CC.PixelIdxList{c}=pixelIndexList{i};
        c=c+1;
    end
end
CC.NumObjects=length(CC.PixelIdxList); %+1  
L=labelmatrix(CC);