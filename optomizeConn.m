function [bw, loc]= optomizeConn_10(gray, L, bias)

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

%%%%%%%%%%%%%%%%%%%%%for testing
% % gray=cir_index;
% gray=outputImage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hist metrics
loc_init= (graythresh(gray)-0.02)*255; % initial guess for optimal thresh
h=histogram(gray(gray>0), 'BinWidth', 1);
% stdv=std(h.Values(1:2:end));
[pks, locs]=findpeaks(smooth(h.Values, 21), 'MinPeakHeight', 1000,...
    'SortStr', 'descend', 'MinPeakProminence', 100);
% pks=pks(1:2); locs=locs(1:2);
pks=pks([1, end]);
locs=locs([1, end]);
% maybe need to use heights, not prom to sort peaks
locs=sort(h.BinEdges(locs));

% select limits and go
% a=round(mean([loc_init, loc_init, locs(1)]));
% b=round(mean([loc_init, loc_init, locs(2)])); 
a=round(0.25*(locs(2)-locs(1))+locs(1));
b=round(-0.25*(locs(2)-locs(1))+locs(2));
connSlope =2; % initialize
bw=gray>0; % initialize
disp('calculating gray con matrix...')
glcms=graycomatrix(gray, 'Offset', [0 1; -1 0],'NumLevels', 256);
disp('Done.')
% imagesc(imadjust((mat2gray(glcms(:,:,1))))); pause(.2)
% figure;

%%
clear Conn per ar level level_prev cc spc MasterMetric
% level=double(max(max(gray)));
level=b+10;
level_prev=level;
c= 0; %counter
while connSlope > 1
    c=c+1;
    level(c)=level_prev-1-1*round((level_prev>b)*(level_prev-b)/4);
    level_prev=level(c);
    bw_prev(:,:,1)=bw;
    bw_prev(:,:,2)=bw_prev(:,:,1); %oldest previous (2 apart)
    bw=gray>level(c);
    imagesc(bw); title(['Level= ', num2str(level(c))]); pause(.01)
    G(c)=getframe(gcf);
    Conn(c)=sum(sum(glcms(level(c):end,level(c):end,1)))...
        +sum(sum(glcms(level(c):end,level(c):end,2)));
%     Conn(c)=sum(sum(bw));
    per(c)=sum(sum(bwperim(bw)));
    ar(c)=sum(sum(bw));
    cc(c)=bwconncomp(bw);
    spc(c)=length(unique(L(bw)));  % # of superpixels classified as water
    disp(['    ', num2str(Conn(c))])
    if (bw_prev(:,:,2)==bw & bw_prev(:,:,1)==bw)
        disp('ended connectivity search bc subsequent images were same')
        break
    elseif level(c) <a % safety strap
        disp('ended connectivity search bc outside search region')
        break
    end
end
% [shldr, loc]=rightShoulder_5(level, Conn, bias); % bias of 2 seems to fix haze prob..

MasterMetric=spc./[cc.NumObjects].*ar./per; % rais something to a power?
[foo, idx]=max(MasterMetric);
loc=level(idx);
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
    level, rescale(Conn./ar),...
    level, rescale([cc.NumObjects]),...
    level, rescale(spc./[cc.NumObjects]),...
    level, rescale(MasterMetric))
hold on; plot(loc, 1.0, 'gV'); hold off
legend({'Area', 'Perim', 'Ar/Per', 'GLCMConn/Ar','ConnComp', 'SPWater/ConnCom', 'SPW/CC*A/P'}, 'Location', 'best')
xlabel('DN threshold'); ylabel('Normalized value');
title('Connectivity Metrics')

