% script to count number of valid observation above grid tiles in dataset
% outputs lists of pairs and also list of unique obs to use for summary
% stats

clear

%% params
useMediumQuality=1;
atUCLA=1;
%% load
if ~isunix
    if atUCLA
        tbl_in='J:\Final\logs\WC_LOG_Summ.xlsx';
    else
        tbl_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\logs\WC_LOG_Summ.xlsx';
    end
else
    tbl_in='/Volumes/Galadriel/Final/logs/WC_LOG_Summ.xlsx';
end

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
if useMediumQuality
    badTilesMask=qualityBad; % change to >2 if want to be more restrictive
else
    badTilesMask=qualityBad | qualityMedium;
end
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
            if useMediumQuality
                if quality(m)==0 || quality(n)==0
                    continue
                end
            else
                if quality(m)==0 || quality(m)==0.5 || quality(n)==0 || quality(n)==0.5
                    continue
                end    
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
if ~isunix
    pairs_path='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\pairs.mat';
else
    pairs_path='/Volumes/Galadriel/Final/analysis/pairs.mat';
    
end
% save(pairs_path, 'pairs');

%% create struct with all centroids and calc stats 
    % for hopes of crunching stats on unpaired lakes
centers=pairs;
repeatIdxs=[];
for i=1:length(pairs)
    repeatIdxs=[repeatIdxs; pairs(i).id'];
end
repeatIdxs=unique(repeatIdxs);
pairedTilesMask=false(size(badTilesMask)); % init
pairedTilesMask(repeatIdxs)=1;
goodUnpairedTilesIdx=find(goodTileMask & ~pairedTilesMask);
goodUnpairedTilesIdx2=setdiff(tiles(goodUnpairedTilesIdx), tiles(pairedTilesMask))

% goodUnpairedTiles=unique(tiles(~pairedTilesMask));
% goodPairedTiles=unique(tiles(pairedTilesMask));

%% create list of tiles to use (2nd try)
% final result: use==1 means use that tile in complete list.  use==2 means
% use an intersection of files.
% total_list.use=zeros(length(dates),1); %init
repeatMask=false(length(dates),1); repeatMask(repeatIdxs)=1;
clear total_list % just in case testing
for i=1:length(dates)
    total_list(i).pair={};
    total_list(i).name=files{i};
    total_list(i).aidx=i;
    if badTilesMask(i) % keep use at 0 if it's a bad tile
        total_list(i).use=0;
        continue
    else
        if ~repeatMask(i) % if good quality and it doesn't repeat, keep it
            total_list(i).use=1;
        else % if good quality and it does repeat, use union of acquisitions (will include geoloc error)
            total_list(i).use=2;
            name=files{i};
            pair_idx=min(find(strcmp({pairs.a}, name) | strcmp({pairs.b}, name)));
            if sum(strcmp(horzcat(total_list.pair), name))==0 % if this tile hasn't already been listed as a pair, use it
                total_list(i).pair=pairs(pair_idx).all';
%                 warning('hit')
            else
                total_list(i).use=3;
            end
        end
    end
end

    % export list
if isunix
    list_out='/Volumes/Galadriel/Final/analysis/total_list.csv';
else
end
total_list_csv=struct2table(total_list);
    % reorder
total_list_csv.zpair=total_list_csv.pair;
total_list_csv.pair=[];

writetable(total_list_csv, list_out)
    % manual removals from pairs list:
