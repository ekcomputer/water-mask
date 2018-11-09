function gorman_thresh = GorminThreshold(hist_counts, wp, df, dyn_range )

% function to find plateau pts of decaying exponential histogram, as
% described in:
% O’Gorman, L. Binarization and Multithresholding of 
% document Images Using Connectivity. CVGIP Graph. Model. Image Process. 
% 56, 494–506 (1994).

% hist_counts is input histogram/PMF (can be offset horizontally)
% plateaus is 1 x n vector of locations (usually just one) to binarize
% image
% wp is sliding window size, expressed as percentage of
% maximum image intensity.  Should be large to reduce noise, but not larger
% than 'min' intensity diff bw levels.
% df is deltaF, or expected flatness deviation as percent of max eul.
% Doesn't change output unless I'm searching for peaks, rather than max.
% Lower values make peaks more distinct, higher values combine peaks.
% dyn_range is aprox dynamic range of image, (close to 256 for uint8) used
% to compute width of sliding window.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%testing

% hist_counts=eul;
% w=11;
% sigma=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global f
    %compute flatness deviation metric
% profl=abs(movsum(hist_counts, w, 'Endpoints', 'fill')-hist_counts);
w=double(wp*0.01*dyn_range/f.cConn); % convert width percentage to width, based on dynamic range, and controlled for step size not equal to 1
d=zeros(w, length(hist_counts)); % init deviation param as matrix
for i=1:length(d) % each level
    for j=max(1,(i-floor(w/2))):min(length(d), (i+floor(w/2))) % sum each window
        d(j,i)=abs(hist_counts(i)-hist_counts(j));
    end
end
d=(sum(d)./sum(d>0)*w)'; % normalize for edges, where window is less than w wide

    %compute flatness deviation metric
sigma=df*0.01*max(d);
profl=exp(-0.5*(d./sigma).^2);

    % compute thresh using simple min
% [~,gorman_thresh]=max(profl);
[~,gorman_thresh]=min(d);
% [~,pks,~,prom]=findpeaks(profl);

    % if using multiple thresh levels:
% thresh_ranges=find(profl>max(profl)*(1-0.01*df))

%% plot test
if f.plot
    figure
    level=1:length(hist_counts); % offset level for plotting purposes
    plot(level, hist_counts/max(hist_counts)*max(d), level, d)
    yyaxis right
    plot(level, profl)
    legend({'eul (norm.)','deviation', 'gaussian profile'}, 'Location', 'southwest')
    ttext{1}=sprintf('O''Gorman Threshold:');
    ttext{2}=sprintf('  wp = %9.1f%%, dF = %9.1f%%', wp, df);
    title(ttext, 'Interpreter', 'tex')
end