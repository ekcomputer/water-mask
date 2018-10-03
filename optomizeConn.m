function [bw, loc]= optomizeConn(gray, L, NoValues, bias)
% V19 fixing issue below again, this time adapting for images with little
% water
% V18 fixes wrong valley issue by sorting by height, not prom!
% V17 works on scaling for black images, etc.
% V16 cleans up mem usage
% V15 inclues masking of image and ignoring high and low hist values
%V 14 includes better mask
% V13 inclues automatic selection of second-highest MasterMetric peak, if
% it is withiin 2% of highest peak and is to the left (lower DN)
% V12 includes sorting by prominence
% V11 removes uncecessary metrics
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
% hist metrics

% put NoData values in gray on the low end
disp('Optimizing Connectivity...')
gray(NoValues)=0;
if 1==1    
    loc_init= (graythresh(gray(~NoValues))-0.02)*255; % initial guess for optimal thresh
    h=histogram(gray(gray>min(gray(:))& gray<max(gray(:))), 'BinWidth', 1, 'BinLimits', [0, max(gray(:))]);
    fig=get(gcf); %record orig figure
    % stdv=std(h.Values(1:2:end));
    [pks, locs, prom]=findpeaks([0; smooth(h.Values, 7); 0], 'MinPeakHeight', 50,...
        'SortStr', 'descend', 'MinPeakProminence', 100,...
        'MinPeakDistance', 100);
    locs=locs-1;
    [pks, heightIndx]=sort(pks.^2.*prom, 'descend'); % or just peaks...
%     if length(pks)==1 & locs(1)<0.5*max(gray(:)) % fixes error if second peak is saturated at right end
%         locs(2)=max(gray(:));
%     if length(pks)==1 & locs(1)>=0.5*max(gray(:))
%         locs(2)=min(gray(:));
%     else
%         pks=pks(heightIndx(1:2));locs=locs(heightIndx(1:2));
%     end
    % maybe need to use heights, not prom to sort peaks
    pks=pks(heightIndx(1:2));locs=locs(heightIndx(1:2));
    locs=sort(h.BinEdges(locs));

    % select limits and go
    % a=round(mean([loc_init, loc_init, locs(1)]));
    % b=round(mean([loc_init, loc_init, locs(2)])); 
    a=min(255, round(locs(1)));
    b=max(0, round(locs(2)));

    if b-a <=4
        a=min(gray(:)); b=max(gray(:));
    end
    connSlope =2; % initialize
    bw=gray>0; % initialize
    % disp('calculating gray con matrix...')
    % glcms=graycomatrix(gray, 'Offset', [0 1; -1 0],'NumLevels', 256);
    % disp('Done.')
    % imagesc(imadjust((mat2gray(glcms(:,:,1))))); pause(.2)
    % figure;

    %%
    clear Conn per ar level level_prev cc spc MasterMetric
    % level=double(max(max(gray)));
    level=b; %10
    level_prev=level;
    c= 0; %counter
    figure;
    while connSlope > 1
        c=c+1;
        level(c)=level_prev-2; %-1*round((level_prev>b)*(level_prev-b)/4);
        level_prev=level(c);
        bw_prev(:,:,1)=bw;
        bw_prev(:,:,2)=bw_prev(:,:,1); %oldest previous (2 apart)
        bw=gray>level(c);
    %     imagesc(bw); title(['Level= ', num2str(level(c))]); pause(.01)
    %     G(c)=getframe(gcf);
    %     Conn(c)=sum(sum(glcms(level(c):end,level(c):end,1)))...
    %         +sum(sum(glcms(level(c):end,level(c):end,2)));
    %     Conn(c)=sum(sum(bw));
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

    MasterMetric=spc./[cc.NumObjects].*ar./per; % rais something to a power?
    % [maxConn, maxConnIdx]=max(MasterMetric);
    %
    try
        [pks, locs]=findpeaks(MasterMetric,... % for MasterMetric
            'SortStr', 'descend', 'MinPeakProminence', 10);
    catch pause
    end
    if isempty(pks) | length(pks)==1 %this block is new
        [pks, locs]=max(MasterMetric);
        pks=[pks,pks]; locs=[locs;locs];
    end
    if length(pks)>=2 % this block is new
        pks=pks(1:2); locs=locs(1:2);
    else 
        pks=[pks,pks]; locs=[locs;locs];
    end
    if (pks(1)-pks(2))/pks(1) <=.02 & level(locs(2))<level(locs(1))
        loc=level(locs(2));
        fprintf('\tNote: second highest connectivity peak automatically chosen\n')
    else
        loc= level(locs(1));
    end
    %%
    % figure
    % plot(level, spc./[cc.NumObjects])
    % title(['Connectivity.  Bias=', num2str(bias)]);
    % xlabel('Threshold level (DN)')
    % ylabel('Connectivity Index')


    bw=gray>loc;
    figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias),...
        ' Level=', num2str(loc)])

    % Test plotting
    figure; plot(level, rescale(ar), level, rescale(per),...
        level, rescale(ar./per),...
        level, rescale([cc.NumObjects]),...
        level, rescale(spc./[cc.NumObjects]),...
        level, rescale(MasterMetric))
    hold on; plot(loc, 1.0, 'gV'); hold off
    legend({'Area', 'Perim', 'Ar/Per','ConnComp', 'SPWater/ConnCom', 'SPW/CC*A/P'}, 'Location', 'best')
    xlabel('DN threshold'); ylabel('Normalized value');
    title('Connectivity Metrics')
    figure(fig.Number); hold on; plot(loc, max(h.Values(:)), 'gV');
    text(loc, 0.9*max(h.Values(:)), {'Connectivity', 'threshold'}, 'color', 'g')
    plot(loc_init,1.0*max(h.Values(:)), 'bV'); hold off
    text(loc_init,0.7*max(h.Values(:)),...
        {'Otsu', 'threshold'}, 'Color', 'b')
else
    loc=im2uint8(graythresh(gray(~NoValues))); % convert to DN
    bw=gray>=loc;
end