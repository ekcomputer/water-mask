function bw_out=sizeFilter(bw, minSize)
% bw_out=sizeFilter(bw, minSize)
% filters regions in binary image below minSize (inclusive)
% uses 8-connectivity
% output is also binary

Lold=bwlabel(bw, 8);
stats=regionprops(bw, 'Area');lstats=length(stats);
allowableAreaIndexes = ([stats.Area] >= minSize);
keeperIndexes = find(allowableAreaIndexes); 
L = ismember(Lold, keeperIndexes); 
bw_out=L>0; clear bw
% L=bwlabel(bw_out);
stats = regionprops(bw_out, 'Area');
fprintf('Size filter completed.  Removed %u regions.  \nNew total= %u regions.\n',...
    lstats-length(stats), length(stats))