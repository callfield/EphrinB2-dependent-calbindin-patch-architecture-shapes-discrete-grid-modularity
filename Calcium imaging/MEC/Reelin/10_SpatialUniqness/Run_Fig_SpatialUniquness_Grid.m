clear all; close all

addpath(pwd);
addpath("function");


load('../Data.mat');

CDir = pwd;

% "Mean_ErrDist","Mean_zCORRMAP","ErrDist","zCORRMAP"
load(strcat(pwd, "\SpatialUniquness_Data.mat"));

th_tags = {'r02p05','r02p01',  'r03p05','r03p01', ...
           'r04p05','r04p01', 'r05p05','r05p01', ...
           'r06p05','r06p01', 'r08p05','r08p01'};  % labels used as field/file suffixes

%% Run the analysis for each threshold tag
for iTag = 12
    tag_name = th_tags{iTag};
    fprintf('===== tag = %s =====\n', tag_name);

    % -------- GROUP1 (cond=1) --------
    zGROUP1_mean   = [];
    zGROUP1_total  = [];
    zGROUP1_max    = [];
    each_zGROUP1_max = nan(2500, 20);   % adjust if the number of sessions changes
    mean_zGROUP1_max = nan(1, 20);
    name_GROUP1   = strings(1, 20);
    n = 1;

    for s = 1:5        % GROUP1 samples
        for t = 1:3    % trial
            DATA = GROUP1{s, t};
            numGrid = length(find(DATA(:, 3) > 0));
            if numGrid > 9 
            
                % Check that the result struct exists.
                if ~isempty(ErrDist{1, s, t}) 
    
                    MEAN_mat  = ErrDist{1, s, t}(iTag).MEAN;
                    TOTAL_mat = ErrDist{1, s, t}(iTag).TOTAL;
                    MAX_mat   = ErrDist{1, s, t}(iTag).MAX;
    
                    zGROUP1_mean  = [zGROUP1_mean;  MEAN_mat(:)];
                    zGROUP1_total = [zGROUP1_total; TOTAL_mat(:)];
                    zGROUP1_max   = [zGROUP1_max;   MAX_mat(:)];
    
                    % Store one vector per session.
                    n_elem = numel(MAX_mat);
                    each_zGROUP1_max(1:n_elem, n) = MAX_mat(:);
                    mean_zGROUP1_max(1, n) = mean(MAX_mat(:), "omitmissing");
                    name_GROUP1(1, n) = strcat(SampleName{1, s}, " t", num2str(t));
    
                    n = n + 1;
                end
            end
        end
    end

    % -------- GROUP2 (cond=2) --------
    zGROUP2_mean   = [];
    zGROUP2_total  = [];
    zGROUP2_max    = [];
    each_zGROUP2_max = nan(2500, 20);
    mean_zGROUP2_max = nan(1, 20);
    name_GROUP2   = strings(1, 20);
    n = 1;

    for s = 1:7
        for t = 1:3
            if s == 1 && t == 3
            else
                DATA = GROUP2{s, t};
                numGrid = length(find(DATA(:, 3) > 0));
                if numGrid > 9
                    if ~isempty(ErrDist{2, s, t}) 
        
                        MEAN_mat  = ErrDist{2, s, t}(iTag).MEAN;
                        TOTAL_mat = ErrDist{2, s, t}(iTag).TOTAL;
                        MAX_mat   = ErrDist{2, s, t}(iTag).MAX;
        
                        zGROUP2_mean  = [zGROUP2_mean;  MEAN_mat(:)];
                        zGROUP2_total = [zGROUP2_total; TOTAL_mat(:)];
                        zGROUP2_max   = [zGROUP2_max;   MAX_mat(:)];
        
                        n_elem = numel(MAX_mat);
                        each_zGROUP2_max(1:n_elem, n) = MAX_mat(:);
                        mean_zGROUP2_max(1, n) = mean(MAX_mat(:), "omitmissing");
                        name_GROUP2(1, n) = strcat(SampleName{2, s}, " t", num2str(t));
        
                        n = n + 1;
                    end
                end
            end
        end
    end

    % Remove unused columns that are all NaN.
    all_nan_GROUP1 = all(isnan(each_zGROUP1_max), 1);
    each_zGROUP1_max(:, all_nan_GROUP1) = [];
    mean_zGROUP1_max(:, all_nan_GROUP1) = [];
    name_GROUP1(:,  all_nan_GROUP1)     = [];

    all_nan_GROUP2 = all(isnan(each_zGROUP2_max), 1);
    each_zGROUP2_max(:, all_nan_GROUP2) = [];
    mean_zGROUP2_max(:, all_nan_GROUP2) = [];
    name_GROUP2(:,  all_nan_GROUP2)     = [];

    % -------- Excel output (one file per tag) --------
    xlsx_file = fullfile(CDir, sprintf('GridMaxFR_MaxErroDist_%s.xlsx', tag_name));

    writematrix(zGROUP1_max,      xlsx_file, 'Sheet', 'GROUP1_all',   "AutoFitWidth", false);
    writematrix(zGROUP2_max,      xlsx_file, 'Sheet', 'GROUP2_all',  "AutoFitWidth", false);
    writematrix(each_zGROUP1_max, xlsx_file, 'Sheet', 'GROUP1',       "AutoFitWidth", false);
    writematrix(each_zGROUP2_max, xlsx_file, 'Sheet', 'GROUP2',      "AutoFitWidth", false);
    writematrix(name_GROUP1,      xlsx_file, 'Sheet', 'GROUP1 name',  "AutoFitWidth", false);
    writematrix(name_GROUP2,      xlsx_file, 'Sheet', 'GROUP2 name', "AutoFitWidth", false);

    %% ==== Figure generation for this tag ====




    % ---- Max (log) ----
    clf
    DATA1 = log10(zGROUP1_max);
    DATA2 = log10(zGROUP2_max);
    XLABEL = "Max Error Distance (log10(cm))";
    TITLE  = sprintf("Max Error distance (RateMap Corr Analysis) [%s]", tag_name);

    Fn_draw_cdf_2sample(DATA1, DATA2, XLABEL, SAMPLELABELS, TITLE)
    exportgraphics(gcf, fullfile(CDir, sprintf("MaxErrorDist_log_cdf_GROUP1vsGROUP2_%s.pdf", tag_name)));
    clf

    BIN = 0.05;
    Fn_draw_histogram_2sample(DATA1, DATA2, XLABEL, SAMPLELABELS, TITLE, BIN)
    exportgraphics(gcf, fullfile(CDir, sprintf("MaxErrorDist_log_hist_GROUP1vsGROUP2_%s.jpg", tag_name)), 'Resolution', 300);
    clf

    % ---- Max (non-log) ----
    clf
    DATA1 = zGROUP1_max;
    DATA2 = zGROUP2_max;
    XLABEL = "Max Error Distance (cm)";
    TITLE  = sprintf("Max Error distance (RateMap Corr Analysis) [%s]", tag_name);

    Fn_draw_cdf_2sample(DATA1, DATA2, XLABEL, SAMPLELABELS, TITLE)
    exportgraphics(gcf, fullfile(CDir, sprintf("MaxErrorDist_nolog_cdf_GROUP1vsGROUP2_%s.jpg", tag_name)), 'Resolution', 300);
    clf

    BIN = 1;
    Fn_draw_histogram_2sample(DATA1, DATA2, XLABEL, SAMPLELABELS, TITLE, BIN)
    exportgraphics(gcf, fullfile(CDir, sprintf("MaxErrorDist_nolog_hist_GROUP1vsGROUP2_%s.jpg", tag_name)), 'Resolution', 300);
    exportgraphics(gcf, fullfile(CDir, sprintf("MaxErrorDist_nolog_hist_GROUP1vsGROUP2_%s.pdf", tag_name)));
    exportgraphics(gcf, ...
    fullfile(CDir, sprintf("MaxErrorDist_nolog_hist_GROUP1vsGROUP2_%s.emf", tag_name)), ...
    "ContentType","vector");
    clf





end


montage_tags = th_tags(12);
montage_tag_name = montage_tags{1};
nrows = 1;
ncols = 1;

% --- Max: log histogram ---
make_tag_montage( ...
    CDir, montage_tags, ...
    'MaxErrorDist_log_hist_GROUP1vsGROUP2_%s.jpg', ...  % source image name for each tag
    sprintf('MaxErrorDist_log_hist_GROUP1vsGROUP2_%s_montage', montage_tag_name), ...
    nrows, ncols);

% --- Max: non-log CDF ---
make_tag_montage( ...
    CDir, montage_tags, ...
    'MaxErrorDist_nolog_cdf_GROUP1vsGROUP2_%s.jpg', ...
    sprintf('MaxErrorDist_nolog_cdf_GROUP1vsGROUP2_%s_montage', montage_tag_name), ...
    nrows, ncols);

% --- Max: non-log histogram ---
make_tag_montage( ...
    CDir, montage_tags, ...
    'MaxErrorDist_nolog_hist_GROUP1vsGROUP2_%s.jpg', ...
    sprintf('MaxErrorDist_nolog_hist_GROUP1vsGROUP2_%s_montage', montage_tag_name), ...
    nrows, ncols);




close all
