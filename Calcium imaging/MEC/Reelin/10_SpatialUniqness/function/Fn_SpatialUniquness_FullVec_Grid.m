function [Result] = Fn_SpatialUniquness_FullVec_Grid(DIR, caFr, r_list, p_list, th_tags, DV_pos, Grid_Scale)

num_SelectCell   = 10;          % number of select cell set
num_SHUFFLE      = 100;


num_SelectFrame  = caFr * 1000; % number of frames used per frame set

load(strcat(DIR, "\ST_dF_grid_aut_data.mat"), 'Grid_Cells');
AnaCellset       = Grid_Cells(:)';  % keep as a row vector
nGrid = numel(AnaCellset);
if num_SelectCell > nGrid
    num_SelectCell = nGrid;
end

Trk = readmatrix(fullfile(DIR, "ST_PCI_Ca_behav_track.csv"));
vlim = 2;
moveFrame = find(Trk(:, 4) > vlim); % velocity over limit (cm/s)
num_mFrame = numel(moveFrame);
nx = 50; % x bin
ny = 50; % y bin
nPix = nx * ny;
nTh = numel(r_list);
pathParts = split(DIR, "\");
samplename = strcat( ...
    char(strrep(pathParts(end-1), '_','-')), " ", ...
    char(strrep(strrep(pathParts(end), '_',' '), 'OF','')) );

CDir = pwd;


%% 
% Build two-dimensional features from Grid_Scale and DV_pos for each grid cell.
% 
% Grid_Scale(i)
% Randomly sample cell sets and score each set by mean pairwise distance.
% Higher scores indicate broader coverage of Grid_Scale and DV_pos.
% Keep the highest-scoring sets as shuffle cell sets.
%
%
%
%
%


nRandomSamples = 20000000;
% Grid_Scale and DV_pos are vectors for all cells in the same order.
[bestSets_local, ~] = select_gridcell_sets( ...
    Grid_Scale, DV_pos, AnaCellset, ...
    num_SelectCell, num_SHUFFLE, nRandomSamples);

% bestSets_local indexes into grid_idx, so convert it to the original cell IDs.
% Store the original cell IDs for each selected cell set.
bestCellSets = AnaCellset(bestSets_local);  % size: [nKeepSets, nCellsPerSet]
% Update the shuffle count after removing duplicate sets.
num_SHUFFLE = size(bestCellSets,1);

%% Define frame sets
num_RAND = ceil(num_mFrame / num_SelectFrame);

if num_RAND == 1
    RAND = 1;  % frames 1 to num_SelectFrame
else
    tmp_merge  = ceil((num_RAND * num_SelectFrame - num_mFrame) / (num_RAND - 1));
    RAND       = 1 + (0:(num_RAND-1)) * (num_SelectFrame - tmp_merge);
end

%% Allocate summary arrays only; do not keep the full CORRMAP in memory.
mean_ErrDist_all  = nan(num_SHUFFLE, num_RAND, nx, ny, nTh);
total_ErrDist_all = nan(num_SHUFFLE, num_RAND, nx, ny, nTh);
max_ErrDist_all   = nan(num_SHUFFLE, num_RAND, nx, ny, nTh);


tic
for ff = 1:num_RAND

    % Target frames
    % mFRAME =  moveFrame(RAND(ff) : RAND(ff) + num_SelectFrame - 1);
    mFRAME = RAND(ff):min(RAND(ff) + num_SelectFrame - 1, num_mFrame);

    % 50 x 50 rate map (cells x x-bin x y-bin)
    tmpGSrate_map = Fn_rate_map_forCorr_Mframe(DIR, mFRAME);

    
    for SHUFFLE = 1:num_SHUFFLE
        cells_idx = bestCellSets(SHUFFLE, :);   % 1 x nCells
    
        % [nCells x nx x ny]
        R = tmpGSrate_map(cells_idx, :, :);
    
        % Flatten to [nCells x nPix].
        % The second-dimension index 1..nPix corresponds to sub2ind([nx,ny],x,y).
        R2 = reshape(R, [], nPix);   % [] = nCells
    
        % Compute all-pixel vs all-pixel correlation matrices.
        % C_full(i,j) = correlation between pixels i and j (pairwise NaN).
        % P_full(i,j) = p-value for that correlation.
        [C_full, P_full] = corr(R2, 'rows', 'pairwise');   % [nPix x nPix]
    
        % % Keep these only if memory allows.
    %     Corr_all{SHUFFLE} = C_full;
    %     Pval_all{SHUFFLE} = P_full;
    % end
    % 


    % for SHUFFLE = 1:num_SHUFFLE
    %     C_full = Corr_all{SHUFFLE};
    %     P_full = Pval_all{SHUFFLE};
    
        for BaseX = 1:nx
            for BaseY = 1:ny
                base_idx = sub2ind([nx, ny], BaseX, BaseY);
    
                corr_map = reshape(C_full(base_idx, :), nx, ny);
                pval_map = reshape(P_full(base_idx, :), nx, ny);
    
                for k = 1:nTh
                    r_th = r_list(k);
                    p_th = p_list(k);
    
                    I2bi = (corr_map >= r_th) & (pval_map < p_th);
    
                    % bwlabel -> regionprops -> distance -> mean/total/max
                    bwl    = bwlabel(I2bi);
                    biarea = regionprops(bwl, 'Area', 'Centroid');
    
                    if isempty(biarea)
                        mean_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k)  = NaN;
                        total_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k) = NaN;
                        max_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k)   = NaN;
                        continue;
                    end
    
                    area_array   = vertcat(biarea.Area);
                    center_array = vertcat(biarea.Centroid);  % (x,y)
    
                    big = (area_array >= 8);
                    if ~any(big)
                        mean_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k)  = NaN;
                        total_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k) = NaN;
                        max_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k)   = NaN;
                        continue;
                    end
    
                    center_big = center_array(big, :);
     %%% CUTION! (BaseY,BaseX)
                    tmp_err_vec = pdist2( ...
                        [BaseY, BaseX], ...
                        [center_big(:,1), center_big(:,2)]);
    
                    mean_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k)  = mean(tmp_err_vec,  "omitnan");
                    total_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k) = sum(tmp_err_vec,   "omitnan");
                    max_ErrDist_all(SHUFFLE, ff, BaseX, BaseY, k)   = max(tmp_err_vec, [], "omitnan");
                end %k
            end %BaseY
        end %BaseX
    end %SHUFFLE



end % ff
toc




for k = 1:nTh
    % ----- Threshold information -----
    r_th = r_list(k);
    p_th = p_list(k);
    tag  = th_tags{k};  % e.g. 'r02p05'

    % === Average across frame sets (dimension 2: num_RAND) ===
    % Extract the 4D array for this threshold k.
    Mean_ErrDist_k  = mean_ErrDist_all( :, :, :, :, k );   % [SHUFFLE, num_RAND, 50, 50]
    Total_ErrDist_k = total_ErrDist_all(:, :, :, :, k);
    Max_ErrDist_k   = max_ErrDist_all(  :, :, :, :, k);

    % Average across the frame-set dimension; keep shape [SHUFFLE,1,50,50].
    
    Mean_ErrDist_meanframeset  = mean(Mean_ErrDist_k,  2, "omitnan"); % [S,1,50,50]
    Total_ErrDist_meanframeset = mean(Total_ErrDist_k, 2, "omitnan"); % [S,1,50,50]
    Max_ErrDist_meanframeset   = mean(Max_ErrDist_k,   2, "omitnan"); % [S,1,50,50]
    
    % For each cell set, sum the 50 x 50 error map and choose the minimum.
   if num_SHUFFLE > 1
    % ---- Mean ----
    score_mean = inf(num_SHUFFLE,1);
    for s = 1:num_SHUFFLE
        M = squeeze(Mean_ErrDist_meanframeset(s,1,:,:));   % [50,50]
        valid = ~isnan(M);
        if any(valid, 'all')
            score_mean(s) = sum(M(valid));  % sum valid pixels only
        end
    end
    [~, idx_best_mean] = min(score_mean);  % cell set with the smallest total error
    Mean_ErrDist_meanframeset_bestcellset = squeeze( ...
        Mean_ErrDist_meanframeset(idx_best_mean,1,:,:) );  % [50,50]

    % ---- Total ----
    score_total = inf(num_SHUFFLE,1);
    for s = 1:num_SHUFFLE
        T = squeeze(Total_ErrDist_meanframeset(s,1,:,:));  % [50,50]
        valid = ~isnan(T);
        if any(valid, 'all')
            score_total(s) = sum(T(valid));
        end
    end
    [~, idx_best_total] = min(score_total);
    Total_ErrDist_meanframeset_bestcellset = squeeze( ...
        Total_ErrDist_meanframeset(idx_best_total,1,:,:) ); % [50,50]

    % ---- Max ----
    score_max = inf(num_SHUFFLE,1);
    for s = 1:num_SHUFFLE
        X = squeeze(Max_ErrDist_meanframeset(s,1,:,:));    % [50,50]
        valid = ~isnan(X);
        if any(valid, 'all')
            score_max(s) = sum(X(valid));
        end
    end
    [~, idx_best_max] = min(score_max);
    Max_ErrDist_meanframeset_bestcellset = squeeze( ...
        Max_ErrDist_meanframeset(idx_best_max,1,:,:) );    % [50,50]

    else
        % If SHUFFLE = 1, use that 50 x 50 map directly and fix index=1.
        Mean_ErrDist_meanframeset_bestcellset  = squeeze(Mean_ErrDist_meanframeset(1,1,:,:));   % [50,50]
        Total_ErrDist_meanframeset_bestcellset = squeeze(Total_ErrDist_meanframeset(1,1,:,:));
        Max_ErrDist_meanframeset_bestcellset   = squeeze(Max_ErrDist_meanframeset(1,1,:,:));
    
        idx_best_mean  = 1;
        idx_best_total = 1;
        idx_best_max   = 1;
    end



%{
    % % Average across the frame-set dimension.
    % Mean_ErrDist_meanframeset  = squeeze(mean(Mean_ErrDist_k,  2, "omitnan")); % [SHUFFLE,50,50] or [50,50]
    % Total_ErrDist_meanframeset = squeeze(mean(Total_ErrDist_k, 2, "omitnan"));
    % Max_ErrDist_meanframeset   = squeeze(mean(Max_ErrDist_k,   2, "omitnan"));
    % 
    % % === Minimum across cell sets (SHUFFLE) ===
    % if num_SHUFFLE > 1
    %     Mean_ErrDist_meanframeset_bestcellset  = squeeze(min(Mean_ErrDist_meanframeset,  [], 1)); % [50,50]
    %     Total_ErrDist_meanframeset_bestcellset = squeeze(min(Total_ErrDist_meanframeset, [], 1));
    %     Max_ErrDist_meanframeset_bestcellset   = squeeze(min(Max_ErrDist_meanframeset,   [], 1));
    % else
    %     Mean_ErrDist_meanframeset_bestcellset  = Mean_ErrDist_meanframeset;   % [50,50]
    %     Total_ErrDist_meanframeset_bestcellset = Total_ErrDist_meanframeset;
    %     Max_ErrDist_meanframeset_bestcellset   = Max_ErrDist_meanframeset;
    % end
%}
    % Final MEAN/TOTAL/MAX for this threshold.
    MEAN  = Mean_ErrDist_meanframeset_bestcellset;
    TOTAL = Total_ErrDist_meanframeset_bestcellset;
    MAX   = Max_ErrDist_meanframeset_bestcellset;

    % ----- Visualization and output (one file per threshold) -----
    th_str  = sprintf(" (r >= %.2f, p < %.3f)", r_th, p_th);
    % Numeric suffix
    th_str2 = sprintf("_r%.2f_p%.3f", r_th, p_th);
    % Use this instead to include the tag name in file names.
    % th_str2 = "_" + tag;
    % === Max Error ===
    imagesc(MAX);
    xticks([]); yticks([]);
    daspect([100 100 100])
    colormap(jet);
    cb = colorbar;
    ylabel(cb, "Max Error Distance (cm)", 'FontSize', 11);
    title({samplename; ...
        strcat("Max error from Top ", num2str(num_SelectCell), ...
               " FR Grid cell", th_str)}, ...
        'FontSize', 14);

    exportgraphics(gcf, ...
        strcat(CDir, "\MaxError_GridmaxFR\Grid MaxError averagedFrameSet maxFR ", ...
               samplename, th_str2, ".jpg"));
    clf

    
    Result(k).tag   = tag;
    Result(k).MEAN  = MEAN;
    Result(k).TOTAL = TOTAL;
    Result(k).MAX   = MAX;

    Result(k).idx_best_mean  = idx_best_mean;
    Result(k).idx_best_total = idx_best_total;
    Result(k).idx_best_max   = idx_best_max;
    
    % ---- Corresponding cell ID sets from bestCellSets ----
    Result(k).bestCellSet_mean  = bestCellSets(idx_best_mean,  :);  % 1 x num_SelectCell
    Result(k).bestCellSet_total = bestCellSets(idx_best_total, :);
    Result(k).bestCellSet_max   = bestCellSets(idx_best_max,   :);



end

close all;

end

