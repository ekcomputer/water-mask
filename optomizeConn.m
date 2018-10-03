% function bw= optomizeConn_6(gray, bias)

% With dynamic programming
% Revised to measure  connectivity with greycomatrix no. of water regions
% connected
% Gray is grayscale image and should be superpixilated and uint8
% Function decreases binary threshold until connectivity is maximized
% bw= output binary classified image

%%%%%%%%%%%%%%%%%%%%%for testing
% % gray=cir_index;
gray=outputImage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b]=deal(30, 140); % bounds for possible threshold.  Tune for better perf.
level=double(max(max(gray)));
level_prev=level;
connSlope =2; % initialize
bw=gray>0; % initialize
disp('calculating gray con matrix...')
% glcms=graycomatrix(gray, 'Offset', [0 1; -1 0],'NumLevels', 256);
disp('Done.')
% imagesc(imadjust((mat2gray(glcms(:,:,1))))); pause(.2)
% figure;

%%
c= 0; %counter
while connSlope > 1
    c=c+1;
    level(c)=level_prev-1-1*round((level_prev>b)*(level_prev-b)/4);
    level_prev=level(c); 
    bw_prev=bw;
    bw=gray>level(c);
    imagesc(bw); title(['Level= ', num2str(level(c))]); pause(.01)
    G(c)=getframe(gcf);
%     Conn(c)=trace(glcms(a:end,a:end,1))+trace(glcms(a:end,a:end,2));
    Conn(c)=sum(sum(bw));
    per(c)=sum(sum(bwperim(bw)));
    ar(c)=sum(sum(bw));
    cc(c)=bwconncomp(bw);
    disp(['    ', num2str(Conn(c))])
    if bw_prev==bw|c == 80 % safety strap
        break
    end
end
[shldr, loc]=rightShoulder_5(level, Conn, bias); % bias of 2 seems to fix haze prob..
plot(level, Conn)
title(['Connectivity.  Bias=', num2str(bias)]); xlabel('Threshold level (DN)')
ylabel('Connectivity Index')
hold on; plot(loc, max(Conn), 'gV'); hold off
bw=gray>loc;
figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias)])
