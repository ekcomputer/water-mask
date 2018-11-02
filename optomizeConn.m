function [bw, loc]= optomizeConn(gray, L, NoValues, bias)

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
otsu=im2uint8(graythresh(gray(~NoValues))); % convert to DN
global f

if 1==1    % branch to use this algorithm
    figure;fig=get(gcf); %record orig figure
    h=histogram(gray(gray>min(gray(:))& gray<max(gray(:))), 'BinWidth',...
        1, 'BinLimits', [0, max(gray(:))]);
    a=15; b=175; % bounds to test
    level=[a:b]; % vector of levels to test
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
            corner = level(RosinThreshold(eul));
            
        catch corner = 256; % in case RosinThresh fails
            warning('RosinThresh failed.')
        end
%     cornerG = GorminThreshold(eul);
        %plotting
    if f.plot
        figure; plot(level, eul/max(eul), level, eul.*per/max(eul.*per),...
            level, eul.*ar/max(eul.*ar))
        legend({'Euler No.', 'E*Per', 'E*Ar'}, 'Location', 'best')
        title('Connectivity Metrics')
        xlabel('DN threshold'); ylabel('Normalized value');
        hold on; plot(corner, 1.0, 'gV'); hold off
        bw=gray>corner;
        figure; imagesc(bw); axis image;
        title(['Initial Mask.  Bias=', num2str(bias),...
        ' Level=', num2str(corner)])

        figure(fig.Number); hold on; plot(corner, max(h.Values(:)), 'gV');
        text(corner, 0.9*max(h.Values(:)), {'Connectivity', 'threshold'}, 'color', 'g')
        plot(otsu,1.0*max(h.Values(:)), 'bV'); hold off
        text(double(otsu),0.8*max(h.Values(:)),...
            {'Otsu', 'threshold'}, 'Color', 'b')
    end
   loc=corner;
   bw=gray>=loc; bw(NoValues)=0;
else % bypass algorithm and just use otsu thresh!
    bw=gray>=otsu;
    loc=otsu;
end
