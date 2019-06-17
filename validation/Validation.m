% computes confusion matrix given input raster files
% written in a very sloppy way...

% restarted 2/2019 for paper
% Revised 2/5/18 to print stats to file for easy use in comparing
% multiple classifications
% Revised 2/22/19 to use standard notation for map-reference confusion
% matrix


%% Choose input classifcation file
clear
saveTables=false; % save file?
saveOutput=false;
n=1;
% im_dir_in='F:\AboveDCSRasterManagement\CanadaAlbersTranslate\';
im_dir_in='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\Final\';
dir_in='D:\ArcGIS\shapes\Validation\rast\';
dir_out='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\validation\';
outputOut=[dir_out, 'output.mat'];
r_out=[dir_out, 'r_class.csv'];
files=cellstr(ls([dir_in, '*.tif']));
lfiles=length(files); 
disp('files found:')
disp(files)
r.out=[]; % init

%%
for n=1:lfiles
        % Choose input 'truth' file
    file_in=[dir_in, files{n}];
    fprintf('Processing file:\n\t%s\n', files{n})
    [bw, bw_R]=geotiffread(file_in);
    disp('Importing file...')
    
        % remove water bodies smaller than 40m
    bw=sizeFilter(bw, 40);
        % Choose input 'classification test' file
    test_name=['WC',files{n}(2:end)];
    test_in=[im_dir_in, test_name];

    [test, test_R]=geotiffread(test_in);
    test=test==1;
        % clip extents to match
    im={~bw, test};
        %% prep to pad smaller file to align pixels in geographic locationshelp
        % find out which file is most NW (assuming cols start from N and rows
        % start from W
    R=[bw_R; test_R];
        % find which image to pad on each side
    [top,north]=min(vertcat(R(:).YWorldLimits)); northFile=north(1); % file with most N upper extent
        top=max(top);
    [left,west]=max(vertcat(R(:).XWorldLimits)); westFile=west(1); % file with most W left extent
        left=min(left);
    [btm,south]=max(vertcat(R(:).YWorldLimits)); southFile=south(2); % file with most W left extent
        btm=min(btm);
    [right,east]=min(vertcat(R(:).XWorldLimits)); eastFile=east(2); % file with most W left extent
        right=max(right);
           
        % compute pad
    for i=1:length(im)
        [NW_intrins(i,1), NW_intrins(i,2)] = worldToIntrinsic(R(i),left, top); % outside bounds. order: x,y
        [SE_intrins(i,1), SE_intrins(i,2)] = worldToIntrinsic(R(i),right, btm);
        buf(i,1)=NW_intrins(i,2)-R(i).YIntrinsicLimits(1); % top buffer
        buf(i,2)=NW_intrins(i,1)-R(i).XIntrinsicLimits(1); % left buffer
        buf(i,3)=-SE_intrins(i,2)+R(i).YIntrinsicLimits(2); % btm buffer
        buf(i,4)=-SE_intrins(i,1)+R(i).XIntrinsicLimits(2); % right buffer
    end
        % pad images
        buf=round(buf);
    for i=1:length(im)
        im{i}=im{i}(1+buf(i,1):end-buf(i,3), 1+buf(i,2):end-buf(i,4));
    end

    %% Create log file
    log=[dir_out, 'IndivLogs\LOG-', files{n}(1:end-4), '.txt'];
    log_ConMat{n}=[dir_out, 'ConMat{n}-', files{n}(1:end-4), '.csv'];
    log_A=[dir_out, 'AccMat-', files{n}(1:end-4), '.csv'];
    % diary([dir_out, 'LOG', files{n}(1:end-4), '.txt'])
    delete(log)
%     diary([dir_out, 'DIARY-', files{n}(1:end-4), '.txt'])
    fid=fopen(log, 'w+');
    fprintf('File:  %s\n', file_in(33:end))
    fprintf(fid, 'File:  %s\n', file_in)
        % update raster values to clips
    bw=im{1}; test=im{2}; clear im
    r.lin=[bw(:), test(:)];
    r.out=[r.out; r.lin];
    diff=double(test)-double(bw);
    figure;
    imagesc(diff); axis image; title({'Difference image', files{n}}); drawnow

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     CF=zeros(3,3); %confusion matrix

    CF(1,1,n)=sum(sum(~test & ~bw));
    CF(1,2,n)=sum(sum(~test & bw));
    CF(2,1,n)=sum(sum(test & ~bw));
    CF(2,2,n)=sum(sum(test & bw));

    CF(3,:,n)=sum(CF(:,:,n),1);
    CF(:,3,n)=sum(CF(:,:,n),2)';

    ConMat{n}=table(CF(:,1,n), CF(:,2,n), CF(:,3,n), 'VariableNames',...
        {'NotWaterValidation','WaterValidation', 'ColumnTotal'},'RowNames',...
        {'NotWaterTest', 'WaterTest', 'RowTotal'})
    if saveTables
        writetable(ConMat{n}, log_ConMat{n});
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=zeros(2,2); %Accuracy{n} matrix
    A(1)=CF(1,1,n)/CF(1,3,n)*100;
    A(2)=CF(2,2,n)/CF(2,3,n)*100;
    A(3)=CF(1,1,n)/CF(3,1,n)*100;
    A(4)=CF(2,2,n)/CF(3,2,n)*100;

    Accuracy{n}=table(A(:,1), A(:,2), 'Rownames',...
        {'NotWater','Water'},'VariableNames',...
        {'ProducerAccuracy','UsersAccuracy'})
    if saveTables
        writetable(Accuracy{n}, log_A);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out(n).UA=A(2,2); % 1-comission error
    out(n).PA=A(2,1); % 1-omission error
    out(n).OA=100*(CF(1,1,n)+CF(2,2,n))/CF(3,3,n); %overall accuracy
    fprintf('Overall Accuracy{n}:  %.2f%%\n', out(n).OA)
    fprintf(fid, 'Overall Accuracy{n}:  %.2f%%\n', out(n).OA);

    out(n).numValPx=sum(sum(test));
    fprintf('Total validation pixels:  %u\n',out(n).numValPx)
    fprintf(fid, 'Total validation pixels:  %u\n',out(n).numValPx);

    out(n).numTestPx=sum(sum(bw));
    fprintf('Total classified pixels:  %u\n\n',out(n).numTestPx)
    fprintf(fid, 'Total classified pixels:  %u\n\n',out(n).numTestPx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    bwstats=regionprops(bw, 'Area', 'Perimeter');
    teststats=regionprops(test, 'Area', 'Perimeter');
    out(n).sizeX=size(bw,1); out(n).sizeY=size(bw,2);
    out(n).numClassWaterBodies=length(bwstats);
    out(n).file=files{n};
    fprintf('Number classified water bodies:  %u\n', out(n).numClassWaterBodies)
    fprintf(fid, 'Number classified water bodies:  %u\n', out(n).numClassWaterBodies);
    
    out(n).numValWaterBodies=length(teststats);
    fprintf('Number validated water bodies:  %u\n', out(n).numValWaterBodies )
    fprintf(fid, 'Number validated water bodies:  %u\n', out(n).numValWaterBodies);

        % perim stats
    out(n).bwP=sum([bwstats.Perimeter]);
    out(n).mapP=sum([teststats.Perimeter]);
    out(n).pDiff= 200*(out(n).bwP-out(n).mapP)/(out(n).bwP+out(n).mapP);
      
        % area stats
    out(n).bwA=sum([bwstats.Area]);
    out(n).mapA=sum([teststats.Area]);
    out(n).aDiff= 200*(out(n).bwA-out(n).mapA)/(out(n).bwA+out(n).mapA);
    
    fprintf('Classified perimeter sum:  %.2u\n', out(n).bwP )
    fprintf(fid, 'Classified perimeter sum:  %.2u\n', out(n).bwP );

    fprintf('Validated perimeter sum:  %.2u\n', out(n).mapP)
    fprintf(fid, 'Validated perimeter sum:  %.2u\n', out(n).mapP);

    fprintf('Percent difference:  %.2f\n', out(n).pDiff);
    fprintf(fid, 'Percent difference:  %.2f\n', out(n).pDiff);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Producgt matrix
    P=zeros(2,2);
    P(1)=CF(7)*CF(3);
    P(2)=CF(8)*CF(3);
    P(3)=CF(7)*CF(6);
    P(4)=CF(8)*CF(6);

    % Kappa Coeff

    ex= trace(P)/sum(sum(P)); %expected chance agreement
    out(n).K=100*(.01*out(n).OA-ex)/(1-ex);
    fprintf('Estimated kappa:  %.2f%%\n', out(n).K);


    fclose('all');
    diary off
end


%% summarize confusion matrix

% for i=1:length(ConMat)
%     CF_all(:,:,i)=table2array(ConMat{i});
% end
% CF_total=sum(CF_all, 3);
CF_total=sum(CF,3);

        % expand matrix
CF_expanded=CF_total;
CF_expanded(4,2)=CF_total(2,2)/CF_total(3,2)*100;
CF_expanded(2,4)=CF_total(2,2)/CF_total(2,3)*100;
CF_expanded(4,1)=CF_total(1,1)/CF_total(3,1)*100;
CF_expanded(1,4)=CF_total(1,1)/CF_total(1,3)*100;
CF_expanded(4,4)=(CF_total(1,1)+CF_total(2,2))/CF_total(3,3)*100;
%% summarize stats

% out(lfiles+1).OA=sum([out.OA].*[out.sizeX].*[out.sizeY])/sum([out.sizeX].*[out.sizeY]);
out(lfiles+1).numValPx=sum([out.numValPx]);
out(lfiles+1).numValWaterBodies=sum([out.numValWaterBodies]);
out(lfiles+1).numTestPx=sum([out.numTestPx]);
out(lfiles+1).numClassWaterBodies=sum([out.numClassWaterBodies]);
out(lfiles+1).file='Summary';
% out(lfiles+1).pDiff=sum(abs([out.pDiff]).*[out.sizeX].*[out.sizeY])/sum([out.sizeX].*[out.sizeY]);
% out(lfiles+1).K=sum([out.K].*[out.sizeX].*[out.sizeY])/sum([out.sizeX].*[out.sizeY]);
% out(lfiles+1).UA=sum([out.UA].*[out.sizeX].*[out.sizeY])/sum([out.sizeX].*[out.sizeY]);
% out(lfiles+1).PA=sum([out.PA].*[out.sizeX].*[out.sizeY])/sum([out.sizeX].*[out.sizeY]);
out(lfiles+1).PA=CF_expanded(4,1);
out(lfiles+1).UA=CF_expanded(1,4);
out(lfiles+1).OA=CF_expanded(4,4);
out(lfiles+1).bwP=sum([out.bwP]);
out(lfiles+1).mapP=sum([out.mapP]);
out(lfiles+1).bwA=sum([out.bwA]);
out(lfiles+1).mapA=sum([out.mapA]);
out(lfiles+1).pDiff=200*(out(lfiles+1).bwP-...
    out(lfiles+1).mapP)/(out(lfiles+1).bwP+out(lfiles+1).mapP);
out(lfiles+1).aDiff=200*(out(lfiles+1).bwA-...
    out(lfiles+1).mapA)/(out(lfiles+1).bwA+out(lfiles+1).mapA);
   % Producgt matrix
    P=zeros(2,2);
    P(1)=CF_total(7)*CF_total(3);
    P(2)=CF_total(8)*CF_total(3);
    P(3)=CF_total(7)*CF_total(6);
    P(4)=CF_total(8)*CF_total(6);
    ex= trace(P)/sum(sum(P)); %expected chance agreement
    out(lfiles+1).K=100*(.01*out(lfiles+1).OA-ex)/(1-ex);

%% save
if saveOutput
    save(outputOut, 'out')
%     writetable(table(r.out(:,1), r.out(:,2)), r_out);
    disp('Output Saved.')
end