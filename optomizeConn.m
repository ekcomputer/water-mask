function bw= optomizeConn_2(gray, bias)

% Gray is grayscale image and should be superpixilated and uint8
% Function decreases binary threshold until connectivity is maximized
% bw= output binary classified image

%%%%%%%%%%%%%%%%%%%%%for testing
% gray=cir_index;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


level=1.0*max(max(gray));
level_prev=level;
c= 0; %counter
connSlope =2; % initialize
bw=gray>0; % initialize
while connSlope > 1
    c=c+1;
    level(c)=0.9*level_prev;
%     level(c)=level_prev-10;
    level_prev=level(c); 
    bw_prev=bw;
    bw=gray>level(c);
    imagesc(bw); title(['Level= ', num2str(level(c))]); pause(.1)
    area=max(sum(bw));
    perim=sum(sum(bwperim(bw)));
    CircPerim=2*pi*sqrt(area/pi) ;
    perSin(c)=perim/CircPerim; % Sinuosity of perim (inverse of former def)
    disp(['    ', num2str(perSin(c))])
    if level(c)<120 & level(c) > 30 % sweet spot for calculating min
        sweetSpot(c)=c;
    end
    if bw_prev==bw|c == 20 % safety strap
        break
    end
end
[ymin, xmin]=u_btm_1(level, perSin, bias) % bias of 2 seems to fix haze prob..
plot(level, perSin)
title(['Connectivity.  Bias=', num2str(bias)]); xlabel('Threshold level (DN)')
ylabel('Connectivity Index')
hold on; plot(xmin, 2*ymin, 'gV'); hold off
bw=gray>xmin;
figure; imagesc(bw); title(['Initial Mask.  Bias=', num2str(bias)])
