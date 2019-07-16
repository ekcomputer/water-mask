% script to add one more observation to waterDistribution.m, based on
% manual clipping of an overlapping region
% TODO: change abun.file to be only filename and not basepath!!
% DEPRICATED - changes added to waterDistribution.m
clear
struct_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\distrib_0.mat';
struct_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\unique\distrib.mat';
minSize=40;%40
load(struct_in);
%%

pths={'D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\straggler\WC_20170806_S01X_Ch045v030_V1_clip.tif'};
i=length(abun)+1; % index of file to append
for j= 1:length(pths)
    disp(i)
    pth_in=pths{j};
    [msk, R]=geotiffread(pth_in);
    NoValues=(msk==-99);
    L=bwlabel(msk);
    init_stats=regionprops(L, 'PixelIdxList', 'Area');
        % calc limnicity (water fraction)
    abun(i).lim=sum(sum(msk==1))/sum(sum(msk>=0));
    abun(i).land=sum(sum(msk==0));
    abun(i).water=sum(sum(msk==1));  
        % buffer all water region
    msk(NoValues)=0;msk=logical(msk);
    SE=strel('diamond',1);
    NoValues_dil=imdilate(NoValues, SE); % msk is now dilated by one
    NoValues_dil(:,1)=1; NoValues_dil(:,end)=1;
    NoValues_dil(1,:)=1; NoValues_dil(end,:)=1;
    
        % remove edge water
    edgeRegions=setdiff(unique(L(NoValues_dil)), 0);
    edgeRegionPx=vertcat(init_stats(edgeRegions).PixelIdxList);
    L(edgeRegionPx)=0;
    
        % size filter
    smallRegions=[init_stats.Area]<minSize;
    smallRegionsPx=vertcat(init_stats(smallRegions).PixelIdxList);
    L(smallRegionsPx)=0;
 
        % calc regionprops
    clear init_stats
    abun(i).stats=regionprops(L>0, 'Area','Perimeter', 'Centroid','MajorAxisLength', 'PixelList');
    if length(abun(i).stats)==0 % if no water in image!
        [abun(i).stats.lat, abun(i).stats.long] = []; % add WGS84 lat/long 
    else
            % plot histogram of stats
        h=histogram([abun(i).stats.Area]/1e6); xlabel('Area ($km^2$)'); ylabel('Count');
        title(pths{j}, 'Interpreter', 'none')
    %         set(gca, 'XScale', 'log', 'YScale', 'log')
        drawnow
        intrinsic=vertcat(abun(i).stats.Centroid);
        clear world
        [world(:,1), world(:,2)]=intrinsicToWorld(R,intrinsic(:,1), intrinsic(:,2));
        gtinfo=geotiffinfo(pth_in);
        mstruct=geotiff2mstruct(gtinfo); % get map projection structure to convert to lat/long
        for k=1:length(abun(i).stats)
            [abun(i).stats(k).lat, abun(i).stats(k).long] = projinv(mstruct,...
                world(k,1), world(k,2)); % add WGS84 lat/long 
            abun(i).stats(k).minBoundRadius=ExactMinBoundCircle(abun(i).stats(k).PixelList);
            abun(i).stats(k).LehnerDevel=(abun(i).stats(k).Area)/(abun(i).stats(k).minBoundRadius^2*pi);
            abun(i).stats(j).SDF=abun(i).stats(j).Perimeter/(2*sqrt(pi*abun(i).stats(j).Area));
        end
        abun(i).freq_min40=length(abun(i).stats)/(abun(i).water+abun(i).land)*100e6;
        abun(i).stats=rmfield(abun(i).stats, 'PixelList');
    end
        % add addiitonal data
    abun(i).file=pths{j};
    abun(i).file_idx=-99;
        % save every ten times
    if mod(i, 10)==0
        save(struct_out, 'abun')
        fprintf('Saving to %s\n', struct_out)
    end
    i=i+1; % in case I'm adding multiple new files
end
save(struct_out, 'abun')
fprintf('Saving to %s\n', struct_out)
disp('Done.')