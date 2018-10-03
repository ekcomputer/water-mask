% function bw= optomizeConn_1(gray)

% Gray is grayscale image and should be superpixilated and uint8
% Function decreases binary threshold until connectivity is maximized

%%%%%%%%%%%%%%%%%%%%%for testing
gray=cir_index;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


level_init=.7*max(max(gray));
level=level_init;

c= 0; %counter
connSlope =2; % initialize
while connSlope > 1
    c=c+1;
    level(c)=.85^c*level_init;
    bw=gray>level(c);
    imagesc(bw); pause(.1)
    area=bwarea(bw);
    perim=bwarea(bwperim(bw));
    CircPerim=2*pi*sqrt(area/pi) ;
    PerSin(c)=perim/CircPerim; % Sinuosity of perim (inverse of former def)
    disp(['    ', num2str(PerSin(c))])
    if c == 20 % safety strap
        break
    end
end
figure; plot(level, PerSin)
[ma, idx]=max(PerSin)
