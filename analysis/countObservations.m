% script to count number of valid observation above grid tiles in dataset

clear
%% load
tbl_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\WC_LOG_Summ.xlsx';
[tbl, tbl_raw]=xlsread(tbl_in, 1);
tbl=tbl(1:330,:);
tbl_raw=(tbl_raw(2:331,:));
idxs=1:330;
dates=str2double(tbl_raw(:, 37));
%% preprocess

leg1Mask=dates<20170800; % july flights
leg2Mask=dates>=20170800; % august flights
    % cloud filter (not used)
hasClouds=tbl(:,23); hasClouds(isnan(hasClouds))=0;
dontRerun=tbl(:,25); dontRerun(isnan(dontRerun))=0;
    % quality filter
quality=tbl(:,22); quality(quality==2)=0.5;
qualityBad=tbl(:,22)==0; 
qualityGood=tbl(:,22)==1;
qualityMedium=tbl(:,22)==2;

badClouds=hasClouds & dontRerun; % flights with clouds that disrupt change deteciton
badTilesMask=qualityBad; % change to >2 if want to be more restrictive
%% compare

tiles=tbl_raw(:,38);
files={tbl_raw{:,2}};
goodTiles=tiles(~badTilesMask);
goodTileIdxs=idxs(~badTilesMask);
goodTileMask=~badTilesMask;
goodTilesLeg1=unique(tiles(leg1Mask & ~badTilesMask));
goodTilesLeg2=unique(tiles(leg2Mask & ~badTilesMask));
goodTileswRepeat=intersect(goodTilesLeg1, goodTilesLeg2);
% goodIdxswRepeat=find(strcmp(tiles,goodTileswRepeat));
goodTileswRepeatMask=false(size(badTilesMask)); % init
for i=1:length(goodTileswRepeat)
    match=find(strcmp(tiles, goodTileswRepeat{i}));
    goodTileswRepeatMask(match)=true;
end
goodTileswRepeatMaskLeg1=goodTileswRepeatMask & leg1Mask;
goodTileswRepeatMaskLeg2=goodTileswRepeatMask & leg2Mask;
%% list files in pair orderings
    % here i is counter for unique tile index, k is counter for unique
    % pairing between file m with tile index i from leg 1 and file n w tile
    % index i from leg 2.  Final pairs struct will contain all possible
    % combs of m and n- some will not overlap, but will be filtered out
    % later.
k=1;
i=1;
while i <= length(goodTileswRepeat)
    goodIdxswRepeat{i}=find(strcmp(tiles,goodTileswRepeat{i}));
    idx_leg1=goodIdxswRepeat{i}(goodTileswRepeatMaskLeg1(goodIdxswRepeat{i}));
    idx_leg2=goodIdxswRepeat{i}(goodTileswRepeatMaskLeg2(goodIdxswRepeat{i}));
    for m = idx_leg1'
        for n= idx_leg2'
            if quality(m)==0 || quality(m)==0
                continue
            end
            pairs(k).idx_leg1=idx_leg1;
            pairs(k).idx_leg2=idx_leg2;
            pairs(k).id=[m, n];
            pairs(k).all=files(goodIdxswRepeat{i})';
            pairs(k).a=char(files{m});
            pairs(k).b=char(files{n});
            pairs(k).quality=min(quality(m), quality(n));
            pairs(k).clouds=max(hasClouds(m), hasClouds(n));
            k=k+1;
        end
    end    
    i=i+1;
end
% sets=sets'

pairs_path='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\pairs.mat';
% save(pairs_path, 'pairs');