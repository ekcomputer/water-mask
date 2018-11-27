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

leg1=dates<20170800; % july flights
leg2=dates>=20170800; % august flights
hasClouds=tbl(:,23); hasClouds(isnan(hasClouds))=0;
dontRerun=tbl(:,25); dontRerun(isnan(dontRerun))=0;
badClouds=hasClouds & dontRerun; % flights with clouds that disrupt change deteciton

%% compare

tiles=tbl_raw(:,38);
files={tbl_raw{:,2}};
goodTiles=tiles(~badClouds);
goodTileIdxs=idxs(~badClouds);
goodTilesLeg1=unique(tiles(leg1 & ~badClouds));
goodTilesLeg2=unique(tiles(leg2 & ~badClouds));
goodTileswRepeat=intersect(goodTilesLeg1, goodTilesLeg2);
% goodIdxswRepeat=find(strcmp(tiles,goodTileswRepeat));

%% list files in pair orderings

for i=1:length(goodTileswRepeat)
    goodIdxswRepeat{i}=find(strcmp(tiles,goodTileswRepeat{i}));
    sets{i}=files(goodIdxswRepeat{i})';
    pairs(i).a=char(sets{i}(1));
    pairs(i).b=char(sets{i}(end));
    pairs(i).all=sets{i};
end
% sets=sets'

pairs_path='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\analysis\pairs.mat';
save(pairs_path, 'pairs');