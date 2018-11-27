function [bw, loc]= optomizeConn(gray, ~, NoValues, bias)
% bias moves thresh in the direction of bias (negative values move thresh
% down)
% Revised to measure  connectivity with greycomatrix no. of water regions
% connected
% Gray is grayscale image and should be superpixilated and uint8
% Function decreases binary threshold until connectivity is maximized
% bw= output binary classified image

% tune smoothing value and prominence value and minPeakDistance

%%%%%%%%%%%%%%%%%%%%%for testing
% % gray=cir_index;
% gray=outputImage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath D:\Dropbox\Matlab\DownloadedCode\RosinThreshold
otsu_thresh=im2uint8(graythresh(gray(~NoValues))); % convert to DN
global f

if 1==1    % branch to use this algorithm
    disp('Calculating global threshold by connectivity')
    figure;fig=get(gcf); %record orig figure
    h=histogram(gray(gray>min(gray(:))& gray<max(gray(:))), 'BinWidth',...
        1, 'BinLimits', [0, max(gray(:))]);
%     a=15; b=175; % bounds to test- 
    level=[f.aConn:f.cConn:f.bConn]; % vector of levels to test
    [eul,per,ar]=deal(zeros(length(level),1)); % init
    for i=length(level):-1:1 % ramp down from bounds [a,b]
        disp(level(i))
        bwt=gray>=level(i); % temp binary image
        eul(i) = bweuler(bwt, 8);
%         per(i)=sum(sum(bwperim(bwt)));
%         ar(i)=sum(sum(bwt));
    end
    
    %compute edge/corner/unimodal hist change pt:
    try
        rosin_thresh = level(RosinThreshold(eul));

    catch rosin_thresh = 256; % in case RosinThresh fails
        warning('RosinThresh failed.')
    end
    
    % compute thresh using O'Gorman method.  Compute dynamic range as input
    % to function.  
    dyn_range=quantile(gray(:),0.99)-quantile(gray(:),0.01);
    gormin_thresh = level(GorminThreshold(eul,f.wp,f.df, dyn_range))
    loc=gormin_thresh; loc=loc+bias;
    bw=gray>=loc; bw(NoValues)=0;
        %plotting
    if f.plot
        figure; plot(level, eul/max(eul), level, eul.*per/max(eul.*per),...
            level, eul.*ar/max(eul.*ar))
%         legend({'Euler No.', 'E*Per', 'E*Ar'}, 'Location', 'best')
        title('Connectivity Metrics')
        xlabel('DN threshold'); ylabel('Normalized value');
        hold on; 
%         plot(gormin_thresh, 1.0, 'gV'); hold off

        plot(rosin_thresh, 1, 'gV');
        text(rosin_thresh, 0.9*1,...
            {'Rosin', 'threshold'}, 'color', 'g')
        plot(otsu_thresh,1.0*1, 'bV'); 
        text(double(otsu_thresh),0.8*1,...
            {'Otsu', 'threshold'}, 'Color', 'b')
        plot(gormin_thresh,1.0*1, 'rV'); hold off
        text(double(gormin_thresh),0.75*1,...
            {'O''Gormin', 'threshold'}, 'Color', 'r')
        
        figure; imagesc(bw); axis image;
        title(['Initial Mask.  Bias=', num2str(bias),...
        ' Level=', num2str(gormin_thresh)])
        hold off

        figure(fig.Number); hold on; 
        plot(rosin_thresh, max(h.Values(:)), 'gV');
        text(rosin_thresh, 0.9*max(h.Values(:)),...
            {'Rosin', 'threshold'}, 'color', 'g')
        plot(otsu_thresh,1.0*max(h.Values(:)), 'bV'); 
        text(double(otsu_thresh),0.8*max(h.Values(:)),...
            {'Otsu', 'threshold'}, 'Color', 'b')
        plot(gormin_thresh,1.0*max(h.Values(:)), 'rV'); hold off
        text(double(gormin_thresh),0.75*max(h.Values(:)),...
            {'O''Gormin', 'threshold'}, 'Color', 'r')
        title('Image histogram')
    end
else % bypass algorithm and just use otsu thresh!
    bw=gray>=otsu_thresh;
    loc=otsu_thresh;
    loc=loc+bias;
end