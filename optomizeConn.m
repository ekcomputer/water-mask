% function bw= optomizeConn_4(gray, bias)

% Revised to measure  connectivity with greycomatrix no. of water regions
% connected
% Gray is grayscale image and should be superpixilated and uint8
% Function decreases binary threshold until connectivity is maximized
% bw= output binary classified image

%%%%%%%%%%%%%%%%%%%%%for testing
% % gray=cir_index;
gray=outputImage;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a, b]=deal(30, 120); % bounds for possible threshold.  Tune for better perf.
level=1.0*max(max(gray));
level_prev=level;
c= 0; %counter
connSlope =2; % initialize
bw=gray>0; % initialize
while connSlope > 1
    c=c+1;
    level(c)=0.90*level_prev;
    level_prev=level(c); 
    bw_prev=bw;
    bw=gray>level(c);
    figure(1)
    imagesc(bw); title(['Level= ', num2str(level(c))]); pause(.1)
    figure(2); imagesc(imadjust((mat2gray(glcms(:,:,1)))))
    glcms=graycomatrix(outputImage.*uint8(bw), 'Offset', [0 1; -1 0],'NumLevels', b-1+1);
%     Conn(c)=trace(glcms(a:end,a:end,1))+trace(glcms(a:end,a:end,2));
    figure(1)
    Conn(c)=sum(sum(glcms(level(c):end,level(c):end,1)))...
        +sum(sum(glcms(level(c):end,level(c):end,2)));
    disp(['    ', num2str(Conn(c))])
    if bw_prev==bw|c == 20 % safety strap
        break
    end
end
[shldr, loc]=rightShoulder_2(level, Conn, bias); % bias of 2 seems to fix haze prob..
plot(level, Conn)
title(['Connectivity.  Bias=', num2str(bias)]); xlabel('Threshold level (DN)')
ylabel('Connectivity Index')
hold on; plot(loc, max(Conn), 'gV'); hold off
bw=gray>loc;
figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias)])
