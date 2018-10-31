function [bw, minvalue]= optomizeConn(sp_mean, sp_text, L, NoValues, bias)

if 1==1  % branch to use this algorithm
    [idx,C,sumd] = kmeans([double(sp_mean(sp_mean>1)), sp_text(sp_mean>1)], 2,...
        'Distance', 'sqeuclidean', 'EmptyAction', 'drop',...
        'Start', [55, 4.5; 200,5], 'MaxIter', 200);
    [sp_bw_idx]=find(logical(idx-1));
    sp_bw=false(size(sp_mean));
    sp_bw(sp_mean>1)=logical(idx-1);
    minvalue=min(sp_mean(sp_bw))+1;
    fprintf('Threshold: %d\n', minvalue)
    disp('Converting to raster...')
    bw=SP_plot_raster(sp_bw, L, 'greaterthan', 0, 'noplot')>0;
    
        % visualize
    figure; yyaxis left; plot(sp_mean(sp_bw), sp_text(sp_bw), 'K.',...
        sp_mean(~sp_bw), sp_text(~sp_bw),'b.')
    hold on; 
    plot(C(:,1), C(:,2), 'gX');     ylabel('Texture')
    yyaxis right; histogram(sp_mean(sp_mean>1))
    hold off
    title('Scatter plot with kmeans centroids')
    xlabel('Water index')
    ylabel('Hist counts')

else % bypass algorithm and just use otsu thresh!
    minvalue=im2uint8(graythresh(sp_mean(sp_mean>0))); % convert to DN
    bw=SP_plot_raster(sp_mean, L, 'greaterthan', minvalue, 'noplot')>0;
    
        %visualize
    figure; histogram(sp_mean(sp_mean>1))
    title('Histogram for otsu thresh')
    xlabel('Water index')
    ylabel('Hist counts')
    hold on; plot(minvalue, 0, 'gV'); hold off
end

figure; imagesc(bw); axis image
