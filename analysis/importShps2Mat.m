% for loading and converting GLWD and my database to .mat
clear
global env
buffer_analysis=0; % load all buffered datasets
glwd_fused_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\dcs_fused_glwd.shp';
dcs_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_combined\WC_complete.shp';
glwd_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\GLWD1and2.shp';
yf_pth='F:\AlaskanLakeDatabase\resource_map_doi_10_5065_D6MC8X5R\data\ek_out\AlaskanLakes_YF.shp';

% hl_fused_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\dcs_fused_hydroLakes.shp';
hl_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\HydroLAKES_intersect.shp';
hl_fused_pth=env.hl_fused_pth;
hl_global_pth='F:\HydroLAKES_polys_v10_shp\HydroLAKES_polys_v10_shp\HydroLAKES_polys_v10.shp';
% out_dir='D:\GoogleDrive\Research\Lake distributions\';
out_dir='D:\GoogleDrive\Research\Lake distributions\';

[~,glwd_fused]=shaperead(glwd_fused_pth);
[~,dcs]=shaperead(dcs_pth);
[~,glwd]=shaperead(glwd_pth);
% [~,yf]=shaperead(yf_pth);
[~,hl_fused]=shaperead(hl_fused_pth);
[~,hl]=shaperead(hl_pth);
[~,hl_global]=shaperead(hl_global_pth);

% pre-process
% glwd=rmfield(glwd, {'Shape_Area', 'Shape_Leng'});

%% atach region name to dcs from fused dataset
% for i=1:length([dcs.Area])
%     lookup=find([fused.Centroid_x]==dcs(i).Centroid_x); &...
%         [fused.Centroid_y]==dcs(i).Centroid_y)
% lookup=find([fused.Shape_Area]==dcs(i).Area &...
%         [fused.Shape_Leng]==dcs(i).Perimeter);
%     dcs(i).Region=fused(lookup).Region;
% end

%% for multiple buffer analysis
if buffer_analysis
    hl_fused_no_buf_pth='D:\ArcGIS\FromMatlab\CIRLocalThreshClas\shp_fused_glwd\dcs_fused_hydroLakes.shp';
    [~,hl_fused]=shaperead(hl_fused_no_buf_pth);
    files=cellstr(ls([env.buffer_dir, '\*sum.dbf']));
    figure; hold on
    j=length(files)+1;
       b(j).buf='0';
       b(j).buf_num=0;
       b(j).count=length([hl_fused.Area]);
       b(j).tot_area=sum([hl_fused.Area]);
       b(j).attr.Area=[hl_fused.Area]; %getfield(hl_fused, 'Area');
    for i=1:j
       if i~=j
           attr=xlsread([env.buffer_dir, '\',files{i}]);
            tmp1=textscan(files{i}, '%s', 'Delimiter', '_');
            b(i).buf=(tmp1{1, 1}{2, 1}  );
            b(i).buf_num=str2double(b(i).buf);
           b(i).count=length(attr(:,3));
           b(i).tot_area=sum(attr(:,3))
           plplot_simple(attr(:,3), 40e-6, 1.1)
       else
           plplot_simple([hl_fused.Area], 40e-6, 1.1)
       end
    end

    % plplot_simple([glwd_fused.Area], 40e-6, 1.1)
    hold off
    legend_txt={b.buf}; %legend_txt{end+1}='0';
    legend(legend_txt)
    figure; scatter([b.buf_num], [b.count]); box on
    xlabel('Buffer length (m)'); ylabel('Number of water bodies')
    figure; scatter([b.buf_num],[b.tot_area])
    xlabel('Buffer length (m)'); ylabel('Total area (km)')
end
%% write
hl_global=rmfield(hl_global, {'Country', 'Continent', 'Poly_src','Grand_id', 'Vol_res',...
    'Vol_src', 'Depth_avg', 'Dis_avg','Res_time', 'Slope_100', 'Wshd_area', 'Pour_lat',...
    'Pour_long', 'Vol_total'});
save([out_dir, 'LakeDatabases.mat'], 'glwd_fused','dcs','glwd', 'hl', 'hl_fused', 'hl_global')
% save([out_dir, 'AlaskaLakes\YukonFlats.mat'], 'yf')