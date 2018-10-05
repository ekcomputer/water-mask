function classified_out=OBIA_BP_Fun(struct_in)

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
close all
tic
addpath D:\Dropbox\Matlab\Above\
global f
f.pArea=1; %pixel area in meters
f.minSize=35; %min water region size (inclusive) in meters squared
f.bounds=[0.9 1.1]; % region growing bounds for regionFill
f.windex='NDWI'; %water index to use
f.satPercent=0.005;
f.Tlim=15; %texture index cutoff
    % ^ must be higher than 10 in order to have any real effect
f.indexShrinkLim=1.3; % max cir_index value (mult by global thresh) for erosion operation
    % ^ 1 or less has no erosion based on value, >1 becomes increasingly
    % discerning
f.sz=100; %target SP size 
f.NDWIWaterLimit=0.001; % global cutoff to determine if tile has only water       <                                                                                    -
f.NDWIWaterAmount=150; % number of pixels above cutoff to show tile has water    <                                                                                    -
f.NDWILandLimit=-0.05; % global cutoff to determine if tile has only land        <                                                                                    -
f.NDWILandAmount=100; % number of pixels above cutoff to show tile has land 

try
    cir=struct_in.data; % for blockproc
catch 
    cir=struct_in; clear struct_in; % for development on single region
end
[cir_index, NoValues, f.waterFlag, f.meanWaterIndex]= BP_loadData(cir, f.windex, 'satPercent', f.satPercent); 

% cir= normalizeImage_old3(cir); % rescale image so that min =0 and max =
if sum(cir_index(:))==0 %if NoData tile
    classified_out=false(size(cir_index));
    disp('Skipping classification')
else
    cir_index(NoValues)=NaN(length(cir_index(NoValues)),1);
    histogram(cir_index(cir_index>0)); title([f.windex,' | saturation: ', num2str(f.satPercent)])
    %% Segmentation
    % optimum scale is about total pixels/450.
    numRows = size(cir,1);
    numCols = size(cir,2);
    totalPix=numRows*numCols;
    f.Diamond=strel('diamond',2); %f is foo structure-no need to save
    % B= entropyfilt(cir_index, f.Diamond.Neighborhood); % looks at 13 neighbors
    fprintf('Target SP size: %d\n', f.sz) 
    disp('Calc superpixels...')
    tic
    [L,N] = superpixels(cir_index,round(totalPix/f.sz), 'Compactness', 0.001, 'Method', 'slic0', 'NumIterations', 10 );
    fprintf('Done.  Average Superpixel size = %0.0f\n', totalPix/N)
    toc
    % figure
    % net = boundarymask(L);
    % imshow(imoverlay(imadjust(cir_index),net,'cyan'),'InitialMagnification',67)
    % imagesc(L)

    %% % Set the value of each pixel in the output image to the mean value
    % of the superpixel region (EK).

    %% % Set the value of each pixel in the output image to the mean value
    % of the superpixel region (EK).

    outputImage = zeros(size(cir_index),'like',cir_index);
    outputText=zeros(size(cir_index),'double');
    sp_mean=zeros(N,1,'like',cir_index); %sp_dev=zeros(N,1,'like',cir_index);
    sp_text=zeros(N,1,'like',cir_index);

    idx = label2idx(L);
    numRows = size(cir_index,1);
    numCols = size(cir_index,2);
    for i = 1:N
        cir_index_Idx = idx{i};
        outputImage(cir_index_Idx) = mean(cir_index(cir_index_Idx));
        outputText(cir_index_Idx) = std(double((cir_index(cir_index_Idx))), 'omitnan');
        sp_mean(i)=mean(cir_index(cir_index_Idx)); %mean value of superpixel
        sp_text(i)=std(double(cir_index(cir_index_Idx)), 'omitNaN'); %mean entropy value of superpixel
        sp_size(i)=sum(sum(double(cir_index(cir_index_Idx))), 'omitNaN');
        %     sp_dev(labelVal)=std(cir_index(redIdx), 'omitNaN'); %std dev value of superpixel
    %     outputImage(greenIdx) = mean(cir_index(greenIdx));
    %     outputImage(blueIdx) = mean(cir_index(blueIdx));
    end 
    clear idx
    % outputText=SP_plot_raster(sp_text, L);
    figure
    % imshow(outputImage,'InitialMagnification',67)
    % imagesc(outputImage);
    % imagesc(imoverlay(cir, boundarymask(L), 'yellow')); axis image

    %plotting
    subplot(311); histogram(sp_mean); title('meanNDWI')
    subplot(312); histogram(sp_text); title('texture')
    subplot(313); histogram(sp_size); title('size'); clear sp_size
    figure
    %% Binary threshold

    % level=graythresh(cir_index);
    % bw=imbinarize(outputImage, level);
    % figure; imagesc(bw); axis image; title('Binary Output Image')
    % close all
    bias=-0; % moves target point left
    if f.waterFlag(1)==1 %1==1  % there is water    
        [bw, f.level]=optomizeConn_25(outputImage, L, NoValues, bias);
    elseif f.waterFlag(1)==0 %there is no water    
        bw=false(size(outputImage));
        fprintf('No water in block\n')
        f.level=9999;
    elseif f.waterFlag(1)==2 %there is only water
        bw=true(size(outputImage));
        bw(NoValues)=0;
        fprintf('No land in block\n')
        f.level=-9999;
    else error('error EK') 
    end
    % figure; imagesc(initialMask(sub.a:sub.b, sub.c:sub.d, :)); axis image
    % title(['Initial Mask.  Bias=', num2str(bias)])

    %% Safety strap
    f.percentWater=sum(bw(:))/sum(~NoValues(:));
    if f.level< 60 & sum(sum(cir_index>f.level))/sum(sum(cir_index>0)) < 0.4
        bw=false(size(bw));
        warning('Caught by safety strap')
    elseif f.meanWaterIndex < -0.07 & f.percentWater>0.4
        bw=false(size(bw));
        f.waterFlag(1)=0;
        disp('No water detected (Safety strap).')
    end
    %% Size filter

    bw=sizeFilter(bw, f.minSize/f.pArea); %minSize given up front
    % size=sz/3?

    %% Visualize
    figure;
    % imagesc(imoverlay(cir, boundarymask(bw), 'blue'));
    imagesc(imoverlay(cir, (bw), 'blue'));
    axis image;title(['Initial Mask.  Bias=', num2str(bias)])
    pause(0.01)

    %% Log
    
    logfile='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\logs\log.txt';
    fid=fopen(logfile, 'a');
    fprintf(fid, '%0.3f\t\t%0.3f\n', f.waterFlag(2),f.waterFlag(3));
    fclose(fid);
    %% Region shrinking (based on entropy)

    disp('Global Shrink...')
    % E_idx_im=im2uint8(rescale(outputEntropy.*double(255-outputImage), 0, 1));
    % bw=bw&E_idx_im<Elim;
    % f.Tlim=6;
    E_idx_mask=outputText<f.Tlim | outputImage>f.indexShrinkLim*f.level; % mask is to keep
    f.szbefore=sum(bw(:));
    % bwnew=bw&~E_idx_mask; 
    bweroded=E_idx_mask&bw; clear bw
    imagesc(bweroded); axis image 
    % imagesc(imoverlay(cir, boundarymask(bwnew), 'yellow')); axis image;
    clear E_idx_mask
    f.szafter=sum(bweroded(:));
    fprintf('Removed %d pixels, or ~%3.0f regions.\n', f.szbefore-f.szafter,...
        (f.szbefore-f.szafter)/f.sz)


    %% Region filling
    clear outputText
    [regiond, Lnew]=regionFill2_1_6(L,bweroded,outputImage, sp_mean, cir_index); 
    % Lnew=bweroded; warning('skipping regionfill') % skip region filling for test
    % [regiond, Lnew]=regionFill3(L,bw,outputImage, sp_mean,...
    %     outputEntropy, f.Tlim, f.bounds, cir_index); 

    %% Fill NaN's surrounded by water
    classified_out=imfillNaN(Lnew>0, NoValues);

    %% visualize
    imagesc(imoverlay(cir, boundarymask(classified_out), 'yellow')); axis image
    title({['Water index cutoff: ', num2str(f.indexShrinkLim),...,...
        ' | Texture index cutoff: ', num2str(f.Tlim),...
        ' |  Mean SP size: ', num2str(round(totalPix/N))],...
        ['Index: ', num2str(f.windex), ' |  Min size: ', num2str(f.minSize),...
        ' |  Growing bounds: ', num2str(f.bounds(1)), ' ',num2str(f.bounds(2))]},...
        'FontSize', 13)
    toc


    %% Save
    disp(f)
%     classified_out=classified_out;

    % dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Intermediate\';
    % name_out=[name_in(1:end-4), '_C','.tif'];
    % 
    % log_out=[dir_out, 'LOG_', name_out(1:end-4), '.txt'];
    % fT=struct2table(f);
    % writetable(fT, log_out);
    % fid=fopen(log_out, 'a');
    % fprintf(fid, '\n\nFile: %s\nCreated: %s\n', [dir_out, name_out], datetime);
    % fprintf(fid, ['Texture index cutoff: ', num2str(f.Tlim)]); fprintf(fid, '\n');
    % fprintf(fid, ['Mean SP size: ', num2str(round(totalPix/N))]); fprintf(fid, '\n');
    % fprintf(fid, ['Index: ', num2str(windex)]); fprintf(fid, '\n');
    % fprintf(fid, ['Min size: ', num2str(f.minSize)]); fprintf(fid, '\n');
    % fprintf(fid, ['Growing bounds: ', num2str(f.bounds(1)), ' ',num2str(f.bounds(2))]);
    % fclose(fid);
    % fprintf('\tFile saved: %s\n', log_out);
end

