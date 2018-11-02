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
addpath D:\Dropbox\Matlab\DownloadedCode\RosinThreshold
otsu=im2uint8(graythresh(gray(~NoValues))); % convert to DN

if 1==1    % branch to use this algorithm
    figure;fig=get(gcf); %record orig figure
    h=histogram(gray(gray>min(gray(:))& gray<max(gray(:))), 'BinWidth',...
        1, 'BinLimits', [0, max(gray(:))]);
    a=100; b=255; % bounds to test
    level=[a:b]; % vector of levels to test
    [eul,per,ar]=deal(zeros(length(level),1)); % init
    for i=length(level):-1:1 % ramp down from bounds [a,b]
        bwt=gray>=i; % temp binary image
        eul(i) = bweuler(bwt, 8);
        per(i)=sum(sum(bwperim(bwt)));
        ar(i)=sum(sum(bwt));
    end
    
        %compute edge/corner/unimodal hist change pt:
    corner = RosinThreshold(eul);
        %plotting
        
    figure; plot(level, eul/max(eul), level, eul.*per/max(eul.*per),...
        level, eul.*ar/max(eul.*ar))
    legend({'Euler No.', 'E*Per', 'E*Ar'}, 'Location', 'best')
    title('Connectivity Metrics')
    xlabel('DN threshold'); ylabel('Normalized value');
    hold on; plot(level(corner), 1.0, 'gV'); hold off
    bw=gray>corner;
    figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias),...
    ' Level=', num2str(corner)])

    figure(fig.Number); hold on; plot(corner, max(h.Values(:)), 'gV');
    text(corner, 0.9*max(h.Values(:)), {'Connectivity', 'threshold'}, 'color', 'g')
    plot(otsu,1.0*max(h.Values(:)), 'bV'); hold off
    text(double(otsu),0.8*max(h.Values(:)),...
        {'Otsu', 'threshold'}, 'Color', 'b')
    
%     figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias),...
%     ' Level=', num2str(loc)])
% 
%     figure; plot(level, ar/max(ar), level, per/max(per),...
%         level, (ar./per)/max(ar./per),...
%         level, [cc.NumObjects]/max([cc.NumObjects]),...
%         level, (spc./[cc.NumObjects])/max(spc./[cc.NumObjects]),...
%         level, MasterMetric/max(MasterMetric), ...
%         level, MM2/max(MM2))
%     hold on; plot(loc, 1.0, 'gV'); hold off
%     legend({'Area', 'Perim', 'Ar/Per','ConnComp', 'SPWater/ConnCom',...
%         'SPW/CC.5*(A/P)2', 'MM2'}, 'Location', 'best')
%     xlabel('DN threshold'); ylabel('Normalized value');
%     title('Connectivity Metrics')
%     
%     figure(fig.Number); hold on; plot(double(loc), max(h.Values(:)), 'gV');
%     text(double(loc), 0.9*max(h.Values(:)), {'Connectivity', 'threshold'}, 'color', 'g')
%     plot(loc_init,1.0*max(h.Values(:)), 'bV'); hold off
%     text(loc_init,0.7*max(h.Values(:)),...
%         {'Otsu', 'threshold'}, 'Color', 'b')
    
        
        
else % bypass algorithm and just use otsu thresh!
    bw=gray>=otsu;
    loc=otus;
end
