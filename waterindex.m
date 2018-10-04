function [cir_index, bw]=waterindex(cir, waterIndex, NoValues)
% cir is 3-band int image
% cir_index is single precision
% bw is binary threshold with certain sensitivity
fprintf('Creating water index: %s\n', waterIndex)
cir=im2single(cir);
switch waterIndex
    case 'NDWI'
        cir_index=(cir(:,:,3)-cir(:,:,1))./(cir(:,:,3)+cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level-.02);
    case 'BNDWI'
%         woij
    case 'MNDWI'
%         cir_ndwi=(cir(:,:,3)-cir(:,:,1))./(cir(:,:,3)+cir(:,:,1));
    case 'NDVI'
%         cir_ndwi=(cir(:,:,2)-cir(:,:,1))./(cir(:,:,2)+cir(:,:,1));
        cir_index=(.5*(cir(:,:,2)+cir(:,:,3))-cir(:,:,1))./(.5*(cir(:,:,2)+cir(:,:,3))+cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level-.02);    
    case 'WI'
        cir_index=(cir(:,:,3)-cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'MAXNDWI'
        cir_index=(max(cir(:,:,2), cir(:,:,3))-cir(:,:,1))./(max(cir(:,:,2), cir(:,:,3))+cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'SPAN'
        u1=[.353, .522, .463]; %Yukon color
        u2=[.106, .192, .212]; %Dark color
        cir_index=-min(spangle_3(u2, cir),spangle_3(u2, cir));
        level=graythresh(cir_index);
        bw=imbinarize(cir_index, level);
    case 'MAXWI'
        cir_index=(max(cir(:,:,2), cir(:,:,3))-cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'MINNDWI'
        cir_index=(min(cir(:,:,2), cir(:,:,3))-cir(:,:,1))./(min(cir(:,:,2), cir(:,:,3))+cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'MINWI'
        cir_index=(min(cir(:,:,2), cir(:,:,3))-cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'BNDWI_3_100'
        cir_index=(100*cir(:,:,3)-cir(:,:,1))./(100*cir(:,:,3)+cir(:,:,1));
%         cir_ndwi=im2single((cir_ndwi+abs(min(min(cir_ndwi))))./(max(max(cir_ndwi))...
%             +abs(min(min(cir_ndwi)))));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'BNDWI_1_100'
        cir_index=(cir(:,:,3)-100*cir(:,:,1))./(cir(:,:,3)+100*cir(:,:,1));
%         cir_ndwi=im2double((cir_ndwi+abs(min(min(cir_ndwi))))./(max(max(cir_ndwi))...
%             +abs(min(min(cir_ndwi)))));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'BNDWI_1_10'
        cir_index=(cir(:,:,3)-10*cir(:,:,1))./(cir(:,:,3)+10*cir(:,:,1));
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'AWEI'
        cir_index=1*cir(:,:,3)-2.5*cir(:,:,2)-3*cir(:,:,1);
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);
    case 'IR'
        cir_index=cir(:,:,1);
        level=graythresh(cir_index(~NoValues));
        bw=imbinarize(cir_index, level);

end