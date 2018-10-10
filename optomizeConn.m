function [bw, loc]= optomizeConn(gray, L, NoValues, bias)
% V27 uses A/p^2 in masterMetric

% works but bad search region detection
% includes rescale function from r2018a
% includes auto optimal region detection and improved metrics
% L is superpixel label matrix
% updated to use better conn metrics
% With dynamic programming
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

%% First 'findpeaks': determine optimal range (using divider)
disp('Optimizing Connectivity...')
    % put NoData values in gray on the low end
gray(NoValues)=0;

if 1==1    % branch to use this algorithm
    loc_init= (graythresh(gray(~NoValues))-0.02)*255; % initial guess for optimal thresh
    h=histogram(gray(gray>min(gray(:))& gray<max(gray(:))), 'BinWidth', 1, 'BinLimits', [0, max(gray(:))]);
    curve=smooth(log(h.Values+1), 5);
    fig=get(gcf); %record orig figure
    figure; plot(curve);
    [pks, locs, prom]=findpeaks([0; curve; 0], 'MinPeakHeight', 0,...
        'SortStr', 'descend', 'MinPeakProminence', 0.01,...
        'MinPeakDistance', 15);
    if length(locs)==1
       bw=false(size(gray));
       loc=-99;
       warning('OptomizeConn: Only one peak found in histogram...')
       return
    end
    locs=locs-1
    divider=137; % DN to split histogram
    l.a= max(locs(locs<=divider));
    g.locsAboveDivider=locs>divider;
    [~,g.PromIndex]=sort(prom);
    [g.MaxPromAboveDivider, g.MaxPromAboveDivLoc]= max(prom(g.locsAboveDivider));
    g.locIndex=[1:length(locs)];
    g.locIndexAboveDivider=g.locIndex([g.locsAboveDivider]);
    l.b=g.locIndexAboveDivider(g.MaxPromAboveDivLoc);
    l.b=locs(g.locIndexAboveDivider(g.MaxPromAboveDivLoc));
    if ~isempty(l.a) & ~isempty(l.b) %condition for peaks on either side of divider
        locs=[l.a, l.b]; 
    else
        locs=locs([g.PromIndex]);
        locs=sort(locs([end,end-1]));
        warning('No peaks found on at least one side of hist divider')
    end
    clear l

    a=locs(1); b=locs(2);
    if a>=b % condition to cover my ass re: previous lines
        a=b-8;
    end

    if b-a <=4
        a=min(gray(:)); b=max(gray(:));
        warning('b-a less than 4')
    end
    
    %% Second findPeaks: 
    clear Conn per ar level level_prev cc spc MasterMetric
    Continue = true ; % initialize
    bw=gray>0; % initialize
    level=b; %10
    level_prev=level;
    c= 0; %counter
    figure;
    while Continue
        c=c+1;
        level(c)=level_prev-6; %-1*round((level_prev>b)*(level_prev-b)/4);
        level_prev=level(c);
        bw_prev(:,:,1)=bw;
        bw_prev(:,:,2)=bw_prev(:,:,1); %oldest previous (2 apart)
        bw=gray>level(c);
        
                % plotting, if desired
    %     imagesc(bw); title(['Level= ', num2str(level(c))]); pause(.01)
    
                % continue
        per(c)=sum(sum(bwperim(bw)));
        ar(c)=sum(sum(bw));
        cc_temp=bwconncomp(bw);
        cc(c).NumObjects=cc_temp.NumObjects;
        spc(c)=length(unique(L(bw)));  % # of superpixels classified as water
        disp(['    ', num2str(level(c))])
        if (bw_prev(:,:,2)==bw & bw_prev(:,:,1)==bw)
            fprintf('ended connectivity search bc subsequent images were same. Bounds: %d %d\n',a,b)
            break
        elseif level(c) <a % safety strap
            fprintf('ended connectivity search bc outside search region.  Bounds: %d %d\n',a,b)
            break
        end
    end
    % [shldr, loc]=rightShoulder_5(level, Conn, bias); % bias of 2 seems to fix haze prob..

    %% compute final answer
    MasterMetric=(spc./[cc.NumObjects]).^(0.5).*(ar./per).^2; % rais something to a power?
    [maxmetric, maxloc]   =   max(MasterMetric);
    loc=level(maxloc);  
    bw=gray>loc;

    %% Plot results!
    figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias),...
        ' Level=', num2str(loc)])
    figure; plot(level, rescale(ar), level, rescale(per),...
        level, rescale(ar./per),...
        level, rescale([cc.NumObjects]),...
        level, rescale(spc./[cc.NumObjects]),...
        level, rescale(MasterMetric))
    hold on; plot(loc, 1.0, 'gV'); hold off
    legend({'Area', 'Perim', 'Ar/Per','ConnComp', 'SPWater/ConnCom', 'SPW/CC.5*(A/P)2'}, 'Location', 'best')
    xlabel('DN threshold'); ylabel('Normalized value');
    title('Connectivity Metrics')
    figure(fig.Number); hold on; plot(double(loc), max(h.Values(:)), 'gV');
    text(double(loc), 0.9*max(h.Values(:)), {'Connectivity', 'threshold'}, 'color', 'g')
    plot(loc_init,1.0*max(h.Values(:)), 'bV'); hold off
    text(loc_init,0.7*max(h.Values(:)),...
        {'Otsu', 'threshold'}, 'Color', 'b')

else % bypass algorithm and just use otsu thresh!
    loc=im2uint8(graythresh(gray(~NoValues))); % convert to DN
    bw=gray>=loc;
end