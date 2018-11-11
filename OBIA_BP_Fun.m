function classified_out=OBIA_BP_Fun(struct_in, log_dir, varargin)

%V8 applies NaN right after SP.  improves NaN fill operatoni at end.
%v7- water land tile test happens 1st
% V6 combines global and local into one function!
%V3 creates empty raster if water flag==0
% OBIA_BP_Fun_2 V2 is for including NaN values

%V 7.6 is debugging memory issues by using older mergeRegSimple (V3.2)
% V7.5 deletes unnces vars to save mem
% V7.4 fixes attempts fixes left bar problem and performs shrink before regionfill
% V7_3 incudes merging
% V7_2 more tuning
% V7_1 - debugging V8
% see V6 for kmeans color stuff
% includes entropy filter
% Rewriting to include region growing/shrinking
% Rewritten to detect SP on masked image
% Next: gradient to detect shoulders
% function labeled=super_1.m(cir, waterIndex)
% script to use OBIA and superpixels to classify open water extent
% TODO: clear bw
%%%%%%%%%%%%%%%%%%%%non- function params
% clear
% close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pause(3)
if isempty(varargin)
    varargin{1}='NULL'; varargin{2}='DATE';
end
close all
ttile=tic;
addpath D:\Dropbox\Matlab\Above\

global f
f.plot=0;
f.pArea=1; %pixel area in meters
f.minSize=40; %min water region size (inclusive) in meters squared
f.bounds=[1.5 2.5]; % region growing bounds for regionFill (coeff for std dev) - the higher, the more it grows
f.windex='NDWI'; %water index to use
f.satPercent= 0.002; %how much to enhance image after initial water index band math
f.Tlim=5.3; %texture index cutoff
    % ^ lower Tlim to erode more heavily (but also remove innner lake pixels)
f.indexShrinkLim=1.5; % max cir_index value (mult by global thresh) for erosion operation
    % ^ 1 or less has no erosion based on value, >1 becomes increasingly
    % discerning
f.sz=100; %target SP size 
% f.NDWIWaterLimit=-0.01; % global cutoff to determine if tile has only water       <                                                                                    -
f.NDWIWaterAmount=0.04; % value of pixels above cutoff to show tile has water    <                                                                                    -
% f.NDWILandLimit=0.01; % global cutoff to determine if tile has only land        <                                                                                    -
f.NDWILandAmount=-0.06; % value of pixels above cutoff to show tile has land 
f.useSafetyStrap=0; %1 to incorporate automated check for bad classification based on % classified
f.minGrowSz=5; % min number of SPs in region to allow regiongrowing (prevents shadow growing)
f.wp=10; %wp is sliding window size for O'Gormin threshold, expressed as percentage
f.df=20; % df is deltaF, or expected flatness deviation as percent of max eul for O'Gormin.  (Doesn't matter for now).
f.aConn=45; % min threshold for O'gormin/Connectivity binarizer
f.bConn=220; % max threshold for O'gormin/Connectivity binarizer
f.cConn=5; % step size for optConn.m
f.growMax=20; % max number of region growing iterations (prevents endless loop)
f.maxStd=0.999; % for region growing: max percential of std-devs for std-dev based growing bounds
f.minAreaFact=300; % number of times to multiply min SP size (in meters) to determine medium high and medium low bounds for initial water determination (higher includes more extreme px)
f.regionsLim= 800; % max number of regions to allow growing.  If larger, assume bad classificatoin and don't grow in order to save time.
try
    cir=struct_in.data; % for blockproc
catch 
    cir=struct_in; clear struct_in; % for developer mode of single tile
    struct_in.blockSize=[-99 -99]; 
    struct_in.imageSize=[-99 -99 -99];
    struct_in.location=[-99 -99];
end
    % pre=populate values for log file
[f.percentWater, f.medWaterWaterIndex, f.meanWaterWaterIndex,...
    f.safetyStrap, f.pregrow, f.postgrow, f.level, f.szbefore,...
    f.szafter, f.origLevel] =deal(-99);

if range(cir(:))==0 %if NoData tile
    classified_out=false(size(cir(:,:,1)));
    NoValues=true(size(cir(:,:,1)));
    f.waterFlag(1)=3;
    disp('NoData tile.')
else
        % load data
    [cir_index, NoValues, f.waterFlag, f.medWaterIndex]= BP_loadData(cir, f.windex, 'satPercent', f.satPercent); 
end

if f.waterFlag(1)==0 %there is no water    
    classified_out=false(size(NoValues));
    fprintf('No water in block\n')
    f.level=-99;
    f.szbefore=-99;f.szafter=-99;
elseif f.waterFlag(1)==2 %there is only water
    bw=true(size(NoValues));
    bw(NoValues)=0;
    fprintf('No land in block\n')
    f.level=-99; 
          % Fill NaN's surrounded by water
    classified_out=imfillNaN_2(bw, NoValues); 
    f.szbefore=-99;f.szafter=-99;
elseif f.waterFlag(1)==3;
    disp('..Skipping classification.')
else
    cir_index(NoValues)=NaN(length(cir_index(NoValues)),1);
%     histogram(cir_index(cir_index>0)); title([f.windex,' | saturation: ', num2str(f.satPercent)])
    %% Segmentation
    % optimum scale is about total pixels/450.
    numRows = size(cir,1);
    numCols = size(cir,2);
    totalPix=numRows*numCols;
    f.Diamond=strel('diamond',2); %f is foo structure-no need to save
    % B= entropyfilt(cir_index, f.Diamond.Neighborhood); % looks at 13 neighbors
    fprintf('Target SP size: %d\n', f.sz) 
    disp('Calc superpixels...')
    tsuperpixels=tic;
    [L,N] = superpixels(cir_index,round(totalPix/f.sz), 'Compactness', 0.001, 'Method', 'slic0', 'NumIterations', 10 );
    L=uint32(L);
    %     L(NoValues)=NaN;
    fprintf('Done.  Average Superpixel size = %0.0f\n', totalPix/N)
    toc(tsuperpixels)
    % figure
    % net = boundarymask(L);
    % imshow(imoverlay(imadjust(cir_index),net,'cyan'),'InitialMagnification',67)
    % imagesc(L)

    %% % Set the value of each pixel in the output image to the mean value
    % of the superpixel region (EK).

    %% % Set the value of each pixel in the output image to the mean value
    % of the superpixel region (EK).
    disp('Converting to vector of superpixel indexes...')
    outputImage = zeros(size(cir_index),'like',cir_index);
%     outputSize = zeros(size(cir_index), 'like', cir_index);
%     outputText=zeros(size(cir_index));
    idx = label2idx(L);
%     N=length(idx);
%     sp_mean=zeros(N,1,'like',cir_index); %sp_dev=zeros(N,1,'like',cir_index);
%     sp_text=zeros(N,1);
%     sp_size=zeros(N,1,'like',cir_index);

%     numRows = size(cir_index,1);
%     numCols = size(cir_index,2);
%     tvect=tic;
    cellmean=@(x)mean(cir_index(x));
    sp_mean=uint8(cellfun(cellmean, idx, 'UniformOutput', true));
    cellentropy=@(x)entropy(cir_index(x));
    sp_text=cellfun(cellentropy, idx, 'UniformOutput', true);
    for i = 1:N
%         cir_index_Idx = idx{i};
%         sp_mean(i)=mean(cir_index(cir_index_Idx)); %mean value of superpixel
        outputImage(idx{i}) = sp_mean(i);  %mean entropy value of superpixel 
%         sp_text(i)=entropy(cir_index(cir_index_Idx)); % how to treat NaN?
    end 
%     toc(tvect)
%     outputImage=SP_plot_raster(sp_mean, L, 'complete');
    clear idx
%     toc(tvect)
    % outputText=SP_plot_raster(sp_text, L);
    % imshow(outputImage,'InitialMagnification',67)
    % imagesc(outputImage);
    % imagesc(imoverlay(cir, boundarymask(L), 'yellow')); axis image

    %plotting
    if f.plot
        figure
        subplot(311); histogram(sp_mean(sp_mean>0)); title('Mean NDWI')
        subplot(312); histogram(sp_text(sp_text>0)); title('NDWI Texture')
    end
    %% Binary threshold

    % level=graythresh(cir_index);
    % bw=imbinarize(outputImage, level);
    % figure; imagesc(bw); axis image; title('Binary Output Image')
    % close all
    bias=-0; % moves target point left
    if f.waterFlag(1)==1 %1==1  % there is water    
        [bw, f.level]=optomizeConn(outputImage, L, NoValues, bias);
    end
    % figure; imagesc(initialMask(sub.a:sub.b, sub.c:sub.d, :)); axis image
    % title(['Initial Mask.  Bias=', num2str(bias)])

    %% Safety strap
    f.percentWater=sum(bw(:))/sum(~NoValues(:));
    f.medWaterWaterIndex=median(cir_index(bw));
    f.meanWaterWaterIndex=mean(cir_index(bw));
    f.safetyStrap=0;
    if f.level< 60 && sum(sum(cir_index>f.level))/sum(sum(cir_index>0)) < 0.4
        warning('Caught by safety strap 1')
        f.safetyStrap=1;
        if f.useSafetyStrap 
            bw=false(size(bw));
            f.waterFlag(1)=0;
            disp('Setting BW to zeros.')
        end
    elseif (f.medWaterIndex < -0.38 ) && f.percentWater>0.3
        warning('No water detected (Safety strap 2).')
        f.safetyStrap=2;
        if f.useSafetyStrap
            bw=false(size(bw));
            f.waterFlag(1)=0;
            f.safetyStrap=2;
            disp('Setting BW to zeros.') 
        end
    elseif f.waterFlag(1)~=2 && (( f.level < 87) && f.percentWater>0.35)
        warning('No water detected (Safety strap 3).')
        f.safetyStrap=3;
        if f.useSafetyStrap
            bw=false(size(bw)); 
            f.waterFlag(1)=0;
            disp('Setting BW to zeros.') 
        end
    end
    %% Size filter
    disp('Size Filter #1')
    bw=sizeFilter(bw, f.minSize/f.pArea); %minSize given up front
    % size=sz/3?

    %% Visualize
    if f.plot
        figure;
        % imagesc(imoverlay(cir, boundarymask(bw), 'blue'));
        imagesc(imoverlay(cir, (bw), 'blue'));
        axis image;title(['Initial Mask.  Bias=', num2str(bias)])
        pause(0.01)
    end
    
 
    %% Condition for local threshold/region growing
    f.nRegions=bweuler(bw);
    if f.nRegions>f.regionsLim
        f.safetyStrap=4; % this is probably a bed thresh and I don't want to grow to save time
        fprintf('Detected %d regions, which is more than %d, so skipping growing to save time.\n', f.nRegions, f.regionsLim)
    elseif f.nRegions>f.regionsLim/2 % condition for maybe bad thresh; reduce grow iterations just to be safe
        f.growMax=round(f.growMax/4);
        fprintf('Detected %d regions, which is more than %d/2, so reducing grow iterations to save time.\n', f.nRegions, f.regionsLim)
    end
    if strcmp(varargin{1}, 'local') && f.waterFlag(1)==1  && sum(bw(:))~=numel(bw) && f.safetyStrap ~=4
 
        %% Region shrinking (based on entropy)
        f.pregrow=sum(bw(:));
        disp('Global Erosion...')
        % E_idx_im=im2uint8(rescale(outputEntropy.*double(255-outputImage), 0, 1));
        % bw=bw&E_idx_im<Elim;
        % f.Tlim=6;
%         E_idx_mask=outputText<f.Tlim & outputImage>f.level;
            % removes SPs with high randomness that are within bw.
            % However, keeps pixels f.indexShrinklim times above the bw
            % level, regardless of randomness.
%         E_idx_mask=outputText<f.Tlim & outputImage>f.level | ...
%             outputImage>f.level*f.indexShrinkLim;
        E_idx_mask=sp_text<f.Tlim & sp_mean>f.level | ...
            sp_mean>f.level*f.indexShrinkLim;
        E_idx_mask=SP_plot_raster(E_idx_mask, L, 'noplot')>0;
        f.szbefore=sum(bw(:));
        % bwnew=bw&~E_idx_mask; 
        bweroded=E_idx_mask&bw; clear E_idx_mask bw
%         imagesc(bweroded); axis image 
        % imagesc(imoverlay(cir, boundarymask(bwnew), 'yellow')); axis image;
        f.szafter=sum(bweroded(:));
        fprintf('Removed %d pixels, or ~%3.0f regions.\n', f.szbefore-f.szafter,...
            (f.szbefore-f.szafter)/f.sz)


        %% Region filling
        clear outputText
        [~, L]=regionFill(L,bweroded,outputImage, sp_mean, sp_text, cir_index); 
        % Lnew=bweroded; warning('skipping regionfill') % skip region filling for test
        % [regiond, Lnew]=regionFill3(L,bw,outputImage, sp_mean,...
        %     outputEntropy, f.Tlim, f.bounds, cir_index); 

        %% Re-apply nodata mask in case SP alg included these regions as water
        L(NoValues)=0;

        %% Size filter
        disp('Size Filter #2')
        classified_out=sizeFilter(L>0, f.minSize/f.pArea); %minSize given up front
        f.postgrow=sum(classified_out(:));
        %% Fill NaN's surrounded by water
        classified_out=imfillNaN(classified_out, NoValues);

        
    else
         % Re-apply nodata mask in case SP alg included these regions as water
        bw(NoValues)=0;
          % Fill NaN's surrounded by water
        classified_out=imfillNaN(bw, NoValues); 
        f.szbefore=-99;f.szafter=-99;
    end
    
    %% visualize
    if f.plot
        figure; imagesc(outputImage); axis image;
        title('Mean water index -segmented')

        figure
        imagesc(imoverlay(cir, boundarymask(classified_out==1), 'yellow')); axis image
        title({['Water index cutoff: ', num2str(f.indexShrinkLim),...,...
            ' | Texture index cutoff: ', num2str(f.Tlim),...
            ' |  Mean SP size: ', num2str(round(totalPix/N))],...
            ['Index: ', num2str(f.windex), ' |  Min size: ', num2str(f.minSize),...
            ' |  Growing bounds: ', num2str(f.bounds(1)), ' ',num2str(f.bounds(2))]},...
            'FontSize', 13)
    end
%     toc
    
end

    % Add mask designation
classified_out=int8(classified_out); % changed these TWO LINES
classified_out(NoValues & ~classified_out)= -99;
%% Save
disp(f)
try disp(struct_in)
catch
end
disp('')
%     classified_out=classified_out;

dir_out=log_dir;
% name_out=[name_in(1:end-4), '_C','.tif'];
% 
name_out=varargin{2};
datecode=varargin{3};
%     log_out=[dir_out, 'LOG_', varargin{2}, char(date), '.csv'];
%     log_out_verbose=[dir_out, 'LOG_X_',...
%         varargin{2}(1:min(length(varargin{2}),30)), '_',char(date), '.csv'];
log_out_verbose = [dir_out, 'LOG_',...
    varargin{3}, '.csv']; 
%     fT=struct2table(f);
%     writetable(fT, log_out);
try fid=fopen(log_out_verbose, 'a');
catch disp('Couldn''t open log');
end
filetext=importdata(log_out_verbose, '\n');
if size(filetext,1)<2
    fprintf(fid, 'File: %s\nCreated: %s\n', [dir_out, name_out], datetime);
    fprintf(fid, 'Filename,PixelArea,MinSize,GrowBound_L,GrowBound_U,WaterIndex,SatPercent,TextureLimit,IndexShrinkLimit,SP_TargetSize,NDWIWaterAmound,NDWILandAmount,WaterFlag,MedianLowIndexes,MedianHighIndexes,MedianIndex,GlobalLevel,PercentWater,MedianWaterIndex,MeanWaterIndex,SizeBeforeShrink,SizeAfterShrink,SafetyStrap,BlockSizeY,BlockSizeX,ImageSizeY,ImageSizeX,ImageSizeZ,LocationY,LocationX,TileMinutes,OrigLevel,SzBfrGrow,SzAftrGrow,wp,df\n');   
end
disp('Tile finished.'); disp(datetime)
elapsedTime=toc(ttile); fprintf('Elapsed time:\t%3.2f minutes\n', ...
    elapsedTime/60);
try
    fprintf(fid, '%s,%9.0f,%9.0f,%9.2f,%9.2f,%s,%9.4f,%d,%9.2f,%9.0f,%9.3f,%9.3f,%9.0f,%9.3f,%9.3f,%9.3f,%9.0f,%9.3f,%9.1f,%9.1f,%9.0f,%9.0f,%d,%d,%d,%d,%d,%d,%d,%d,%9.2f,%9.2f,%d,%d,%9.1f,%9.1f\n' ,...
    name_out ,f.pArea,f.minSize,f.bounds(1),f.bounds(2),f.windex,f.satPercent,...
    f.Tlim,f.indexShrinkLim,f.sz,f.NDWIWaterAmount,f.NDWILandAmount,...
    f.waterFlag(1), f.waterFlag(3),f.waterFlag(2),f.medWaterIndex,...
    f.level, f.percentWater,f.medWaterWaterIndex,f.meanWaterWaterIndex,...
    f.szbefore, f.szafter, f.safetyStrap,...
    struct_in.blockSize(1),...
    struct_in.blockSize(2), struct_in.imageSize(1), struct_in.imageSize(2),...
    struct_in.imageSize(3), struct_in.location(1), struct_in.location(2),...
    elapsedTime/60, f.origLevel, f.szbefore, f.szafter, f.wp, f.df);
catch

    try fprintf(fid, '%s, -99\n' , name_out); % defensively record at least name of file
    catch disp('NO LOG FILE WRITTEN...')
    end
end
fclose(fid);
fprintf('\tParam Log File saved: %s\n', log_out_verbose);
