% ========================================================================
% Bayesian Position Decoding Pipeline
% ========================================================================
% This script performs Bayesian decoding of the animal's spatial position
% from hippocampal CA1 calcium activity recorded during open-field
% exploration (100 × 100 cm arena).
%    The goal is to quantify spatial coding accuracy and compare decoding
% performance between WildType and IslandKilled animals.
%    The decoder uses a Gaussian likelihood model and evaluates decoding
% accuracy using cross-validation.
%
% ------------------------------------------------------------------------
% OVERVIEW OF THE ANALYSIS PIPELINE
% ------------------------------------------------------------------------
% 1. Load session data
%    - Calcium activity (CaTrace)
%    - Behavioral tracking (x,y position and speed)
%    - Cell classification tables (place cells, reliable cells, all cells)
% 2. Preprocess neural activity
%    - Align behavioral tracking to calcium timestamps
%    - Restrict analysis to movement epochs (speed > 2 cm/s)
%    - Z-score calcium traces for each cell
% 3. Spatial binning
%    - The arena is divided into 2.5 cm spatial bins
%    - Each frame is assigned to a spatial bin based on animal position
% 4. Bayesian decoder (5-fold cross-validation)
%    - Data are split into 5 contiguous temporal folds
%    - For each fold:
%      TRAIN (80% of frames)
%      - Estimate per-cell mean activity in each spatial bin
%      - Estimate variance of activity in each spatial bin
%      - Apply Gaussian spatial smoothing to reduce noise (sigma=2 bins)
%      TEST (20% of frames)
%      - For each frame, compute the log-likelihood of neural activity for every spatial bin using a Gaussian likelihood model:   log P(activity | position)
%      - Combine with occupancy prior:   P(position | activity) ∝ P(activity | position) · P(position)
%      - Decode position as the spatial bin with maximum likelihood
% 5. Decoding accuracy
%    - Compute Euclidean distance between decoded and true position
%    - Report mean decoding error across frames
% 6. Cell population comparisons
%    Decoding is performed for three cell populations:
%      • Place cells  (• All reliable rate-map cells  • All recorded cells)
% 7. Cell-count matched control analysis
%    WildType sessions often contain more cells than IslandKilled sessions.
%    To control for this:
%      - Both sessions are randomly subsampled to a fixed cell count
%      - The decoder is re-run multiple times (default: 50 iterations)
%      - Mean/median decoding error across subsamples is reported
%    This tests whether decoding differences are driven by cell number.
% 8. Shuffle control
%    As a chance-level control:
%      - Neural activity is circularly time-shifted
%      - Decoding is repeated on shuffled data
%      - This estimates chance decoding accuracy
% 9. Statistical comparisons will be conducted on BayesianDecoding_Hisa_260316_CellCountSweep_Chance_forNaoki.m
%    The current script only reports comparison with full-population and specific subsampled population.
% 10. Output
%    The script saves:
%      • Per-session decoding results (without shuffle, including shuffle makes peformance heavy)
%      • Cell-count matched decoding results
%      • Shuffle controls
% ------------------------------------------------------------------------
% DECODER MODEL
% ------------------------------------------------------------------------
% The decoder assumes Gaussian activity distributions per cell and bin:
%   log P(activity | bin) =
%      -0.5 * Σ_cells [ (x - μ)^2 / σ^2 + log(σ^2) ]
% where
%   x  = observed activity
%   μ  = mean activity in that spatial bin
%   σ² = variance of activity in that spatial bin
% The posterior probability is computed as:
%   P(bin | activity) ∝ P(activity | bin) · P(bin)
% where P(bin) is the spatial occupancy prior.
% ========================================================================




%%
clear
close all

%  Load metrics for place cell / reliable cell / all cell
Cinfo = load("H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\stats_260206\ExportedData\TraceData_AllSuccessfulTrials\PlaceCellInfo_260218.mat");

for i = 1:size(Cinfo.CellList_category_all,1)
    for j = 1:size(Cinfo.CellList_category_all,2)
        if isempty(Cinfo.CellList_category_all{i,j}) == 0
            Cinfo.CellList_category_all{i,j}.Cell_All = Cinfo.CellListT_all_animals{i,j};
        end
    end
end


count = 0;
for CellN_Target = [5:1:60] % number of random subsampling cell numbers
    disp(CellN_Target)
    count = count + 1; % full-population decoding is conducted only at 1st loop

    %% parameters for decoding
    close all

    % Load data of three successful trials
    data_directory   = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\stats_260206\ExportedData\TraceData_AllSuccessfulTrials';

    % Set output directory
    results_directory = 'H:\experiments H drive\BayesResults_260317_BorderCenter';
    results_directory = fullfile(results_directory, strcat(num2str(CellN_Target), 'cells'));

    if ~exist(results_directory, 'dir'), mkdir(results_directory); end

    spatial_bin_size_cm     = 2.5; % arena is divided into 2.5 cm spatial bins
    minimum_speed_cm_per_s = 2; % Restrict analysis to movement epochs (speed > 2 cm/s)
    arena_boundaries_cm    = [0 100; 0 100];   % [xmin xmax; ymin ymax]
    smoothing_sigma_bins   = 2; % gaussian smoothing bins

    number_of_cv_folds        = 5;  %(number of cross-validation. 5-fold: 80% training, 20% test)
    number_of_subsample_iters = 50; % random subsampling is re-run multiple times (default: 50 iterations)
    minimum_cell_count        = 5;  % Only sessions contains >5 cells are used.

    nShuffle = 20; % shuffling iterations, defalt: 20

    border_width_cm = 20; % define border [cm]

    population_field_names = {'Cell_reliableRatemap_PlaceCell', 'Cell_reliableRatemap', 'Cell_All'};
    population_display_names = {'PlaceCells', 'AllReliable', 'AllCells'};
    number_of_populations = length(population_field_names);

    random_seed               = 42;


    kernel_halfwidth = ceil(smoothing_sigma_bins * 3) * 2 + 1;
    gaussian_kernel  = fspecial('gaussian', kernel_halfwidth, smoothing_sigma_bins);


    %% Identify sessions and cells

    CellList_category_all = Cinfo.CellList_category_all;
    mice_str = Cinfo.mice_str;

    for row = 1:size(CellList_category_all, 1)
        for col = 1:size(CellList_category_all, 2)
            entry = CellList_category_all{row, col};
            if ~isempty(entry) && isstruct(entry)
                field_list = fieldnames(entry);
                fprintf('  CellList{%d,%d} fields:\n', row, col);
                for fi = 1:length(field_list)
                    this_field = entry.(field_list{fi});
                    if istable(this_field)
                        fprintf('    .%-40s %dx%d table\n', field_list{fi}, height(this_field), width(this_field));
                    end
                end
                for pop = 1:number_of_populations
                    if isfield(entry, population_field_names{pop})
                        this_table = entry.(population_field_names{pop});
                        if istable(this_table) && any(strcmp(this_table.Properties.VariableNames, 'CellInd'))
                            fprintf('  OK: %s has CellInd (%d cells)\n', population_field_names{pop}, height(this_table));
                        end
                    end
                end
                fprintf('\n'); break;
            end
        end
        if ~isempty(entry) && isstruct(entry), break; end
    end


    fprintf('=== Discovering session files ===\n');
    available_files = dir(fullfile(data_directory, 'CaTrace_And_VideoTrack_*.mat'));
    number_of_files = length(available_files);

    session_list = struct('animal_name', {}, 'recording_day', {}, 'filename', {}, ...
        'mouse_index', {}, 'group_label', {},...
        'CellN_All', {}, 'CellN_Reliable', {}, 'CellN_Place', {});

    for file_idx = 1:number_of_files
        this_filename = available_files(file_idx).name;
        parsed_tokens = regexp(this_filename, 'CaTrace_And_VideoTrack_(.+)_d(\d+)\.mat', 'tokens');
        if isempty(parsed_tokens), continue; end

        animal_name   = parsed_tokens{1}{1};
        recording_day = str2double(parsed_tokens{1}{2});
        mouse_index   = find(strcmp(mice_str, animal_name));
        if isempty(mouse_index), continue; end
        if recording_day > size(CellList_category_all, 2), continue; end

        cell_info_entry = CellList_category_all{mouse_index, recording_day};
        if isempty(cell_info_entry) || ~isstruct(cell_info_entry), continue; end

        % Check if population has enough cells
        has_enough_cells = false;
        for pop = 1:number_of_populations
            if isfield(cell_info_entry, population_field_names{pop})
                this_table = cell_info_entry.(population_field_names{pop});
                if istable(this_table) && height(this_table) >= minimum_cell_count
                    has_enough_cells = true; break;
                end
            end
        end
        if ~has_enough_cells, continue; end

        if startsWith(animal_name, 'n'),     group_label = 'WildType';
        elseif startsWith(animal_name, 'p'), group_label = 'IslandKilled';
        else,                                group_label = 'unknown';
        end

        cellNall = height(cell_info_entry.Cell_uncategorized) + height(cell_info_entry.Cell_unstableRatemap) + height(cell_info_entry.Cell_reliableRatemap);
        cellNreliable = height(cell_info_entry.Cell_reliableRatemap);
        cellNplace = height(cell_info_entry.Cell_reliableRatemap_PlaceCell);

        session_list(end+1) = struct('animal_name', animal_name, ...
            'recording_day', recording_day, 'filename', this_filename, ...
            'mouse_index', mouse_index, 'group_label', group_label, ...
            'CellN_All', cellNall, 'CellN_Reliable', cellNreliable, 'CellN_Place', cellNplace); %#ok<SAGROW>
    end

    number_of_sessions      = length(session_list);
    number_of_wildtype      = sum(strcmp({session_list.group_label}, 'WildType'));
    number_of_islandkilled  = sum(strcmp({session_list.group_label}, 'IslandKilled'));
    fprintf('  Found %d sessions (%d WildType, %d IslandKilled)\n\n', ...
        number_of_sessions, number_of_wildtype, number_of_islandkilled);



    %%  RUN DECODER ON EACH SESSION (full-population)
    %  For each session: Load calcium traces and behavioral tracking and  Z-score traces, restrict to movement epochs
    %     Run 5-fold cross-validated Bayesian decoding , plots areSaved as scatter
    % full-population decoding is conducted only at 1st loop, pre-processed data is temporaly saved as cached data
    if count == 1
        [all_decoding_results, cached_session_data] = func_BayesDec(number_of_populations, number_of_sessions, results_directory, session_list, ...
            data_directory, minimum_speed_cm_per_s, arena_boundaries_cm, spatial_bin_size_cm, gaussian_kernel, number_of_cv_folds, ...
            minimum_cell_count, CellList_category_all, population_field_names, population_display_names, border_width_cm);
        all_decoding_results_temp = all_decoding_results;
    elseif count >1
        all_decoding_results =  all_decoding_results_temp;
    end



    %% 3. MATCHED CELL-COUNT random SUBSAMPLING
    fprintf('\n=== Matched cell-count subsampling ===\n');
    subsampling_timer = tic;

    for pop = 1:number_of_populations
        this_population_label = population_display_names{pop};
        has_result = ~cellfun(@isempty, all_decoding_results{pop});
        if sum(has_result) == 0, continue; end
        valid_results = safe_cat_structs(all_decoding_results{pop}(has_result));

        is_wildtype     = strcmp({valid_results.group_label}, 'WildType'); % Control
        is_islandkilled = strcmp({valid_results.group_label}, 'IslandKilled'); % Casp

        if sum(is_islandkilled)==0 || sum(is_wildtype)==0
            fprintf('  [%s] Need both groups, skip\n', this_population_label); continue;
        end

        islandkilled_cell_counts = [valid_results(is_islandkilled).number_of_cells];
        target_cell_count = CellN_Target;

        fprintf('  [%s] Subsampling WildType to %d cells (median IslandKilled count)\n', ...
            this_population_label, target_cell_count);

        if target_cell_count < minimum_cell_count
            fprintf('  [%s] Target too small, skip\n', this_population_label); continue;
        end

        sessions_with_results = find(has_result);

        for vi = 1:length(sessions_with_results)

            session_idx = sessions_with_results(vi);
            this_result = all_decoding_results{pop}{session_idx};

            if this_result.number_of_cells <= target_cell_count

                this_result.matched_median_error = NaN;
                this_result.matched_mean_error   = NaN;
                this_result.matched_std_error    = NaN;
                this_result.matched_target_cell_count = target_cell_count;
                this_result.matched_all_errors   = NaN;
                % add here other results

                all_decoding_results{pop}{session_idx} = this_result;
                continue;
            end

            %  Retrieve cached session data
            this_cache = cached_session_data{session_idx};

            if isempty(this_cache)
                fprintf('  [%s] %s day%d: cache miss\n', ...
                    this_population_label, this_result.animal_name, this_result.recording_day);
                continue;
            end

            cached_calcium_traces     = this_cache.calcium_traces_raw;
            cached_movement_indices   = this_cache.movement_frame_indices;
            cached_x_bins             = this_cache.x_bin_during_movement;
            cached_y_bins             = this_cache.y_bin_during_movement;
            cached_true_x             = this_cache.true_position_x;
            cached_true_y             = this_cache.true_position_y;
            cached_fold_assignment    = this_cache.cv_fold_assignment;
            cached_frame_interval     = this_cache.frame_interval_sec;
            cached_n_x_bins           = this_cache.number_of_x_bins;
            cached_n_y_bins           = this_cache.number_of_y_bins;
            cached_x_centers          = this_cache.x_bin_centers;
            cached_y_centers          = this_cache.y_bin_centers;

            cached_is_border_frame      = this_cache.is_border_frame; %index in movement
            border_idx = find(cached_is_border_frame);
            center_idx = find(~cached_is_border_frame);

            subsample_timer_this = tic;

            subsample_median_errors = zeros(number_of_subsample_iters, 1);
            subsample_median_errors_border = zeros(number_of_subsample_iters,1);
            subsample_median_errors_center = zeros(number_of_subsample_iters,1);

            original_cell_ids   = this_result.cell_ids;
            original_cell_count = this_result.number_of_cells;

            subsample_per_frame_errors  = cell(number_of_subsample_iters, 1);
            % subsample_decoded_xy        = cell(number_of_subsample_iters, 1);
            subsample_per_frame_errors_border = cell(number_of_subsample_iters, 1);
            subsample_per_frame_errors_center = cell(number_of_subsample_iters, 1);

            subsample_median_errors_shuf = cell(number_of_subsample_iters, 1);
            subsample_median_errors_shuf_med = zeros(number_of_subsample_iters, 1);
            subsample_per_frame_errors_shuf  = cell(number_of_subsample_iters, 1);
            subsample_decoded_xy_shuf        = cell(number_of_subsample_iters, 1);

            subsample_median_errors_shuf_border_med = zeros(number_of_subsample_iters, 1);
            subsample_median_errors_shuf_center_med = zeros(number_of_subsample_iters, 1);

            subsample_per_frame_errors_shuf_border  = cell(number_of_subsample_iters, 1);
            subsample_per_frame_errors_shuf_center  = cell(number_of_subsample_iters, 1);


            % Precompute z-scored traces once per session
            z_all = cached_calcium_traces;

            subset_mean = mean(z_all,1,'omitnan');
            subset_std  = std(z_all,0,1,'omitnan');
            subset_std(subset_std < 1e-10) = 1;

            z_all = (z_all - subset_mean) ./ subset_std;
            z_all(isnan(z_all)) = 0;
            z_all = single(z_all);

            % movement frames pre-extracted
            z_movement = z_all(cached_movement_indices,:);

            % local copies for parfor
            z_all_local = z_all;
            z_movement_local = z_movement;
            movement_idx = cached_movement_indices;
            border_idx_movementInAll = movement_idx(border_idx);
            center_idx_movementInAll = movement_idx(center_idx);
            % ==========================

            rng(random_seed)
            cell_permutations = zeros(number_of_subsample_iters, target_cell_count);
            for iter = 1:number_of_subsample_iters
                cell_permutations(iter,:) = original_cell_ids(randperm(original_cell_count, target_cell_count));
            end

            % parfor iter = 1:number_of_subsample_iters
            for   iter = 1:number_of_subsample_iters
                rng(random_seed + iter) % this did not properly work in parfor loop

                % Randomly pick target_cell_count cells from this session
                random_cell_subset = cell_permutations(iter,:);
                subset_zscored = z_movement_local(:, random_cell_subset);

                % decoding
                [this_median_error, this_mean_error, this_per_frame_errors, this_decoded_x_all, this_decoded_y_all, true_x_out, true_y_out] ...
                    = run_bayesian_decoder(subset_zscored, cached_true_x, cached_true_y, cached_x_bins, cached_y_bins, ...
                    cached_n_x_bins, cached_n_y_bins, cached_x_centers, cached_y_centers, cached_frame_interval, gaussian_kernel, cached_fold_assignment, number_of_cv_folds);

                subsample_median_errors(iter) = this_median_error;
                subsample_per_frame_errors{iter} = single(this_per_frame_errors);


                % ==========================================================
                % ===== Shuffle calc. ======
                % ==========================================================
                shuffle_error = NaN(nShuffle,1);
                shuffle_per_frame_errors = cell(nShuffle,1);
                decoded_xy_shuf = cell(nShuffle,1);


                % parfor i = 1:nShuffle
                for i = 1:nShuffle
                    rng(random_seed + iter*1000 + i);

                    % === All ===
                    nFrames = size(z_all_local,1);
                    shift = randi([round(nFrames*0.1) round(nFrames*0.9)]);
                    % z_shifted = circshift(z_all_local(:, random_cell_subset), shift);
                    z_shifted = z_all_local(:, random_cell_subset);
                    for c = 1:size(z_shifted,2)
                        shift = randi([round(nFrames*0.1) round(nFrames*0.9)]);
                        z_shifted(:,c) = circshift(z_shifted(:,c), shift);
                    end
                    subset_zscored_shuf = z_shifted(movement_idx,:);

                    [median_error, mean_error, per_frame_errors, decoded_x, decoded_y, observed_x, observed_y] ...
                        = run_bayesian_decoder(subset_zscored_shuf, cached_true_x, cached_true_y,cached_x_bins, cached_y_bins, ...
                        cached_n_x_bins, cached_n_y_bins, cached_x_centers, cached_y_centers, cached_frame_interval, gaussian_kernel, cached_fold_assignment, number_of_cv_folds);
                    shuffle_error(i) = median_error;
                    shuffle_per_frame_errors{i} = single(per_frame_errors);

                end

 
                tmp = cat(2, shuffle_per_frame_errors{:});
                shuffle_per_frame_errors = mean(tmp, 2, 'omitnan');
                tmp = cat(2, shuffle_per_frame_errors_border{:});
                shuffle_per_frame_errors_border = mean(tmp, 2, 'omitnan');
                tmp = cat(2, shuffle_per_frame_errors_center{:});
                shuffle_per_frame_errors_center = mean(tmp, 2, 'omitnan');

                subsample_median_errors_shuf_med(iter) = mean(shuffle_error,'omitmissing');
                subsample_median_errors_shuf_border_med(iter) = mean(shuffle_error_border,'omitmissing');
                subsample_median_errors_shuf_center_med(iter) = mean(shuffle_error_center,'omitmissing');
                subsample_per_frame_errors_shuf{iter} = shuffle_per_frame_errors;
                subsample_per_frame_errors_shuf_border{iter} = shuffle_per_frame_errors_border;
                subsample_per_frame_errors_shuf_center{iter} = shuffle_per_frame_errors_center;

                % saving this is too heavy
                % subsample_per_frame_errors_shuf{iter} = per_frame_errors_shuf;
                % subsample_decoded_xy_shuf{iter} = decoded_xy_shuf;

            end

            % subsample_per_frame_errors = cell2mat(subsample_per_frame_errors')';

            elapsed_subsample = toc(subsample_timer_this);

            this_result.matched_median_error = median(subsample_median_errors);
            this_result.matched_mean_error   = mean(subsample_median_errors);
            this_result.matched_std_error    = std(subsample_median_errors);
            this_result.matched_all_errors   = subsample_median_errors;
            this_result.matched_target_cell_count = target_cell_count;
            this_result.matched_per_frame_errors        = cell2mat(subsample_per_frame_errors');
 
            %
            this_result.shuf_matched_median_error = median(subsample_median_errors_shuf_med);
            this_result.shuf_matched_all_errors   = subsample_median_errors_shuf_med;

            this_result.shuf_matched_all_errors_border   = subsample_median_errors_shuf_border_med;
            this_result.shuf_matched_all_errors_center   = subsample_median_errors_shuf_center_med;

            this_result.shuf_matched_per_frame_errors_shuf          = cell2mat(subsample_per_frame_errors_shuf');
            this_result.shuf_matched_per_frame_errors_shuf_border   = cell2mat(subsample_per_frame_errors_shuf_border');
            this_result.shuf_matched_per_frame_errors_shuf_center   = cell2mat(subsample_per_frame_errors_shuf_center');


            all_decoding_results{pop}{session_idx} = this_result;

            fprintf('  [%s] %s day%d: full=%.1f cm, matched(%d cells)=%.1f +/- %.1f cm (%.0f sec)\n', ...
                this_population_label, this_result.animal_name, this_result.recording_day, ...
                this_result.median_decoding_error, target_cell_count, ...
                this_result.matched_median_error, this_result.matched_std_error, elapsed_subsample);

            % add shuffle per-frame mean and std

        end
    end

    fprintf('Subsampling complete in %.0f sec\n', toc(subsampling_timer));


    %%  STATISTICAL COMPARISONS (not used)
    fprintf('\n=== Statistical comparisons ===\n');

    for pop = 1:number_of_populations
        this_population_label = population_display_names{pop};
        fprintf('\n--- %s ---\n', this_population_label);
        has_result = ~cellfun(@isempty, all_decoding_results{pop});
        if sum(has_result)==0, fprintf('  No results\n'); continue; end
        results = safe_cat_structs(all_decoding_results{pop}(has_result));

        is_wildtype     = strcmp({results.group_label}, 'WildType');
        is_islandkilled = strcmp({results.group_label}, 'IslandKilled');

        % Full population comparison
        wildtype_errors     = [results(is_wildtype).median_decoding_error];
        islandkilled_errors = [results(is_islandkilled).median_decoding_error];
        [pvalue_full, ~] = ranksum(wildtype_errors, islandkilled_errors);
        fprintf('  Full population (all cells used):\n');
        fprintf('    WildType:      median=%.1f cm, IQR=%.1f\n', median(wildtype_errors), iqr(wildtype_errors));
        fprintf('    IslandKilled:  median=%.1f cm, IQR=%.1f\n', median(islandkilled_errors), iqr(islandkilled_errors));
        fprintf('    Wilcoxon rank-sum p = %.4f\n', pvalue_full);

        % Matched cell-count comparison
        wildtype_matched     = [results(is_wildtype).matched_median_error];
        islandkilled_matched = [results(is_islandkilled).matched_median_error];
        try
            [pvalue_matched, ~] = ranksum(wildtype_matched, islandkilled_matched);
            fprintf('  Matched cell count:\n');
            fprintf('    WildType:      median=%.1f cm, IQR=%.1f\n', median(wildtype_matched), iqr(wildtype_matched));
            fprintf('    IslandKilled:  median=%.1f cm, IQR=%.1f\n', median(islandkilled_matched), iqr(islandkilled_matched));
            fprintf('    Wilcoxon rank-sum p = %.4f\n', pvalue_matched);
        catch
        end

        % Mixed-effects models
        try
            lme_table = table();
            lme_table.DecodingError = [results.median_decoding_error]';
            lme_table.Group         = categorical({results.group_label})';
            lme_table.Animal        = {results.animal_name}';
            lme_table.CellCount     = [results.number_of_cells]';

            lme_group_only = fitlme(lme_table, 'DecodingError ~ Group + (1|Animal)');
            fprintf('  LME: DecodingError ~ Group + (1|Animal)\n');
            fprintf('    Group effect: beta=%.2f, p=%.4f\n', ...
                lme_group_only.Coefficients.Estimate(2), lme_group_only.Coefficients.pValue(2));

            lme_with_cellcount = fitlme(lme_table, 'DecodingError ~ Group + CellCount + (1|Animal)');
            fprintf('  LME: DecodingError ~ Group + CellCount + (1|Animal)\n');
            fprintf('    Group:     beta=%.2f, p=%.4f\n', ...
                lme_with_cellcount.Coefficients.Estimate(2), lme_with_cellcount.Coefficients.pValue(2));
            fprintf('    CellCount: beta=%.2f, p=%.4f\n', ...
                lme_with_cellcount.Coefficients.Estimate(3), lme_with_cellcount.Coefficients.pValue(3));
        catch error_info
            fprintf('  LME failed: %s\n', error_info.message);
        end
    end

    %% Figures to generate

    wildtype_color     = [0.2 0.4 0.8];
    islandkilled_color = [0.8 0.2 0.2];
    for pop = 1:number_of_populations
        this_population_label = population_display_names{pop};
        has_result = ~cellfun(@isempty, all_decoding_results{pop});

        if sum(has_result)==0, continue; end
        results = safe_cat_structs(all_decoding_results{pop}(has_result));

        is_wildtype     = strcmp({results.group_label}, 'WildType');
        is_islandkilled = strcmp({results.group_label}, 'IslandKilled');
        wildtype_errors          = [results(is_wildtype).median_decoding_error];
        islandkilled_errors      = [results(is_islandkilled).median_decoding_error];
        wildtype_matched_errors  = [results(is_wildtype).matched_median_error];
        islandkilled_matched_errors = [results(is_islandkilled).matched_median_error];
        wildtype_cell_counts     = [results(is_wildtype).number_of_cells];
        islandkilled_cell_counts = [results(is_islandkilled).number_of_cells];

        wildtype_matched_shuf_errors  = [results(is_wildtype).shuf_matched_median_error];
        islandkilled_matched_shuf_errors = [results(is_islandkilled).shuf_matched_median_error];

        summary_figure = figure('Position', [50 50 1800 900], 'Color', 'w');
        sgtitle(sprintf('Bayesian Decoding: %s', this_population_label), 'FontSize', 14, 'FontWeight', 'bold');

        %Full population - all cells used
        subplot(2,4,1);
        group_data = [wildtype_errors(:); islandkilled_errors(:)];
        group_labels = [ones(length(wildtype_errors),1); 2*ones(length(islandkilled_errors),1)];
        boxplot(group_data, group_labels, 'Labels', {'WildType','IslandKilled'}, 'Colors', [0 0 0], 'Symbol', '', 'Widths', 0.5);
        hold on;
        scatter(1+0.15*randn(length(wildtype_errors),1), wildtype_errors, 40, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7);
        scatter(2+0.15*randn(length(islandkilled_errors),1), islandkilled_errors, 40, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);
        ylabel('Median decoding error (cm)');
        [pv,~] = ranksum(wildtype_errors, islandkilled_errors);
        title(sprintf('Full population (p=%.4f)', pv));

        % Matched cell count
        subplot(2,4,2);
        try
            group_data = [wildtype_matched_errors(:); islandkilled_matched_errors(:)];
            group_labels = [ones(length(wildtype_matched_errors),1); 2*ones(length(islandkilled_matched_errors),1)];
            boxplot(group_data, group_labels, 'Labels', {'WildType','IslandKilled'}, 'Colors', [0 0 0], 'Symbol', '', 'Widths', 0.5);
            hold on;
            scatter(1+0.15*randn(length(wildtype_matched_errors),1), wildtype_matched_errors, 40, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7);
            scatter(2+0.15*randn(length(islandkilled_matched_errors),1), islandkilled_matched_errors, 40, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);
            ylabel('Median decoding error (cm)');
            [pv,~] = ranksum(wildtype_matched_errors, islandkilled_matched_errors);
            title(sprintf('Matched cell count (p=%.4f)', pv));
        catch
        end

        % Error vs cell count
        subplot(2,4,3);
        scatter(wildtype_cell_counts, wildtype_errors, 50, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7); hold on;
        scatter(islandkilled_cell_counts, islandkilled_errors, 50, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);
        xlabel('Number of cells'); ylabel('Median decoding error (cm)');
        legend('WildType', 'IslandKilled', 'Location', 'northeast');
        title('Decoding error vs cell count');

        % Full vs matched error (not important)
        subplot(2,4,5);
        try
            scatter(wildtype_errors, wildtype_matched_errors, 50, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7); hold on;
            scatter(islandkilled_errors, islandkilled_matched_errors, 50, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);
            plot([0 60],[0 60],'k--');
            xlabel('Full population error (cm)'); ylabel('Matched cell count error (cm)');
            legend('WildType', 'IslandKilled', 'Location', 'northwest');
            title('Effect of cell count matching'); axis square;
        catch
        end

        %  Cumulative distribution
        subplot(2,4,6);
        [cdf_wt_y, cdf_wt_x] = ecdf(wildtype_errors);
        [cdf_ik_y, cdf_ik_x] = ecdf(islandkilled_errors);
        plot(cdf_wt_x, cdf_wt_y, '-', 'Color', wildtype_color, 'LineWidth', 2); hold on;
        plot(cdf_ik_x, cdf_ik_y, '-', 'Color', islandkilled_color, 'LineWidth', 2);
        xlabel('Median decoding error (cm)'); ylabel('Cumulative fraction');
        legend('WildType', 'IslandKilled', 'Location', 'southeast');
        title('Cumulative distribution of errors');

        % Per-animal summary (each animal = one data, median of median)
        subplot(2,4,4);
        unique_animals = unique({results.animal_name}, 'stable');
        per_animal_median_error = zeros(length(unique_animals), 1);
        per_animal_group        = cell(length(unique_animals), 1);
        for animal_idx = 1:length(unique_animals)
            belongs_to_this_animal = strcmp({results.animal_name}, unique_animals{animal_idx});
            per_animal_median_error(animal_idx) = median([results(belongs_to_this_animal).median_decoding_error]);
            per_animal_group{animal_idx} = results(find(belongs_to_this_animal, 1)).group_label;
        end
        wildtype_animal_medians     = per_animal_median_error(strcmp(per_animal_group, 'WildType'));
        islandkilled_animal_medians = per_animal_median_error(strcmp(per_animal_group, 'IslandKilled'));
        group_data = [wildtype_animal_medians; islandkilled_animal_medians];
        group_labels = [ones(length(wildtype_animal_medians),1); 2*ones(length(islandkilled_animal_medians),1)];
        boxplot(group_data, group_labels, 'Labels', {'WildType','IslandKilled'}, 'Colors', [0 0 0], 'Symbol', '', 'Widths', 0.5);
        hold on;
        scatter(1+0.15*randn(length(wildtype_animal_medians),1), wildtype_animal_medians, 50, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7);
        scatter(2+0.15*randn(length(islandkilled_animal_medians),1), islandkilled_animal_medians, 50, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);
        ylabel('Median decoding error (cm)');
        [pv,~] = ranksum(wildtype_animal_medians, islandkilled_animal_medians);
        title(sprintf('Per-animal median (p=%.4f)', pv));

        % Matched cell count _ shuffled
        subplot(2,4,7);
        try
            %%
            group_data = [wildtype_matched_errors(:); islandkilled_matched_errors(:); wildtype_matched_shuf_errors(:); islandkilled_matched_shuf_errors(:)];
            group_labels = [ones(length(wildtype_matched_errors),1); 2*ones(length(islandkilled_matched_errors),1); ...
                3*ones(length(wildtype_matched_shuf_errors),1); 4*ones(length(islandkilled_matched_shuf_errors),1)];
            boxplot(group_data, group_labels, 'Labels', {'WildType','IslandKilled', 'ContShuf', 'CaspShuf'}, 'Colors', [0 0 0], 'Symbol', '', 'Widths', 0.5);
            hold on;
            scatter(1+0.15*randn(length(wildtype_matched_errors),1), wildtype_matched_errors, 40, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7);
            scatter(2+0.15*randn(length(islandkilled_matched_errors),1), islandkilled_matched_errors, 40, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);
            scatter(3+0.15*randn(length(wildtype_matched_shuf_errors),1), wildtype_matched_shuf_errors, 40, wildtype_color, 'filled', 'MarkerFaceAlpha', 0.7);
            scatter(4+0.15*randn(length(islandkilled_matched_shuf_errors),1), islandkilled_matched_shuf_errors, 40, islandkilled_color, 'filled', 'MarkerFaceAlpha', 0.7);

            ylabel('Median decoding error (cm)');
            [pv,~] = ranksum(wildtype_matched_errors, islandkilled_matched_errors);
            title(sprintf('Matched cell count (p=%.4f)', pv));
        catch
        end


        %%
        saveas(summary_figure, fullfile(results_directory, sprintf('BayesDecoding_%s.png', this_population_label)));
        saveas(summary_figure, fullfile(results_directory, sprintf('BayesDecoding_%s.fig', this_population_label)));
    end

    %% save results


    save(fullfile(results_directory, 'bayesian_decoding_results.mat'), ...
        'all_decoding_results', 'session_list', 'population_display_names', 'population_field_names', ...
        'number_of_subsample_iters', 'number_of_cv_folds', 'spatial_bin_size_cm', 'smoothing_sigma_bins', ...
        'minimum_speed_cm_per_s', 'minimum_cell_count', '-v7.3');

    csv_file = fopen(fullfile(results_directory, 'bayesian_decoding_summary.csv'), 'w');
    fprintf(csv_file, 'Population,Animal,Day,Group,CellCount,MedianError_cm,MeanError_cm,MatchedError_cm,MatchedStdError_cm,MatchedTargetCellCount\n');
    for pop = 1:number_of_populations
        has_result = ~cellfun(@isempty, all_decoding_results{pop});
        valid_indices = find(has_result);
        for j = 1:length(valid_indices)
            r = all_decoding_results{pop}{valid_indices(j)};
            fprintf(csv_file, '%s,%s,%d,%s,%d,%.2f,%.2f,%.2f,%.2f,%d\n', ...
                population_display_names{pop}, r.animal_name, r.recording_day, r.group_label, ...
                r.number_of_cells, r.median_decoding_error, r.mean_decoding_error, ...
                r.matched_median_error, r.matched_std_error, r.matched_target_cell_count);
        end
    end
    fclose(csv_file);

    fprintf('\nDone! All results saved to:\n  %s\n', results_directory);

end



%% functions
%% Bayesian decoder

function [median_error, mean_error, per_frame_errors, decoded_x_all, decoded_y_all, true_x_out, true_y_out] ...
    = run_bayesian_decoder(zscored_traces, true_x, true_y, x_bin_indices, y_bin_indices, number_of_x_bins, number_of_y_bins, ...
    x_bin_centers, y_bin_centers, frame_interval, gaussian_kernel, fold_assignment, number_of_folds)

%  Gaussian likelihood Bayesian decoder with k-fold cross-validation.
% TRAIN: 80% Estimate per-cell mean and variance of z-scored activity in
%           each spatial bin using sparse matrix accumulation.
% SMOOTH: Gaussian-smooth occupancy, mean, and variance maps using
%            imfilter on the full 3D stack , (can  also use conv2 per cell) .
% TEST: 20%  Compute log-likelihood of each spatial bin for all test frames
%           simultaneously via matrix multiplication; decoded position = argmax bin.
%  Returns: median/mean Euclidean error, per-frame errors, decoded positions.

number_of_frames = size(zscored_traces,1);
number_of_cells  = size(zscored_traces,2);
total_spatial_bins = number_of_x_bins * number_of_y_bins;

per_frame_errors = NaN(number_of_frames,1);
decoded_x_all    = NaN(number_of_frames,1);
decoded_y_all    = NaN(number_of_frames,1);
true_x_out       = true_x;
true_y_out       = true_y;

for fold = 1:number_of_folds

    is_test_frame  = fold_assignment == fold;
    is_train_frame = ~is_test_frame;

    nTrain = sum(is_train_frame);
    nTest     = sum(is_test_frame);

    if nTest == 0 || nTrain < 50
        continue
    end

    %  this is to gather activity per spatial bin in training
    training_traces = zscored_traces(is_train_frame,:);
    training_x_bins = x_bin_indices(is_train_frame);
    training_y_bins = y_bin_indices(is_train_frame);

    has_valid_bin = training_x_bins>0 & training_y_bins>0;

    linear_bin_index = sub2ind([number_of_x_bins number_of_y_bins], ...
        training_x_bins(has_valid_bin), training_y_bins(has_valid_bin));

    valid_training_traces = training_traces(has_valid_bin,:);

    %   sparse indiacator matrix
    nValid = size(valid_training_traces,1);
    bin_indicator = sparse( (1:nValid)', linear_bin_index, 1, nValid, total_spatial_bins);

    % occupancy
    occupancy_count = full(sum(bin_indicator,1))';
    % activity sums
    sum_of_activity = bin_indicator' * valid_training_traces;
    % squared activity sums
    sum_of_activity_squared = bin_indicator' * (valid_training_traces.^2);

    occupancy_map       = reshape(occupancy_count, number_of_x_bins, number_of_y_bins);
    mean_activity       = reshape(sum_of_activity, number_of_x_bins, number_of_y_bins, number_of_cells);
    squared_activity    = reshape(sum_of_activity_squared, number_of_x_bins, number_of_y_bins, number_of_cells);

    % ==============================
    % smoothing
    % ==============================
    % Smooth rate maps are created from the traces
    smoothed_occupancy = imfilter(double(occupancy_map), gaussian_kernel,'same','conv');
    smoothed_occupancy(smoothed_occupancy < 1) = NaN;

    mean_activity = imfilter(mean_activity, gaussian_kernel, 'same','conv');
    squared_activity = imfilter(squared_activity, gaussian_kernel, 'same','conv');

    mean_activity = mean_activity ./ smoothed_occupancy;
    squared_activity = squared_activity./smoothed_occupancy;

    variance_map = squared_activity - mean_activity.^2;
    variance_map(variance_map < 0.01) = 0.01;

    mean_per_bin = reshape(mean_activity,total_spatial_bins,number_of_cells);
    variance_per_bin = reshape(variance_map,total_spatial_bins,number_of_cells);

    was_visited = occupancy_count>0;

    % ==============================
    % likelihood terms
    % ==============================
    % Decoding: vectorized log-likelihood via matrix multiply
    %  log P(activity | bin) = -0.5 * sum_cells[(x - mu)^2 / var + log(var)]
    %  I expand the quadratic term
    %    sum[(x - mu)^2 / var] = x^2*(1/var) - 2*x*(mu/var) + mu^2/var
    % the first two terms become matrix multiplies across all test frames all at once.
    %  the third term is constant across frames and precomputed per bin (saving time).

    inverse_variance = 1 ./ variance_per_bin;
    mean_over_variance = mean_per_bin .* inverse_variance;

    per_bin_constant = 0.5 * sum(mean_per_bin.^2 .* inverse_variance,2) ...
        + 0.5 * sum(log(variance_per_bin),2);

    test_frame_traces = zscored_traces(is_test_frame,:);
    test_true_x = true_x(is_test_frame);
    test_true_y = true_y(is_test_frame);

    % Two matrix multiplies: (nTest x nCells) * (nCells x nBins) -> nTest x

    quadratic_term = test_frame_traces.^2 * inverse_variance';
    linear_term    = test_frame_traces * mean_over_variance';

    log_likelihood = -0.5 * quadratic_term + linear_term - per_bin_constant';

    % ==============================
    % PRIOR
    % ==============================
    %P(x∣activity)∝P(activity∣x)⋅P(x)

    smoothed_occ_vector = reshape(smoothed_occupancy,total_spatial_bins,1);
    prior = smoothed_occ_vector ./ sum(smoothed_occ_vector,'omitnan');
    log_prior = log(prior + eps);
    log_likelihood = log_likelihood + log_prior';
    log_likelihood(:,~was_visited) = -Inf;

    [~,most_likely_bin] = max(log_likelihood,[],2);

    % Convert bin indices back to spatial positions (may be prone to small errors, % but not much
    [decoded_x_bin,decoded_y_bin] = ind2sub([number_of_x_bins number_of_y_bins],most_likely_bin);

    decoded_x_position = x_bin_centers(decoded_x_bin)';
    decoded_y_position = y_bin_centers(decoded_y_bin)';

    euclidean_error = sqrt((decoded_x_position(:)-test_true_x(:)).^2 + ...
        (decoded_y_position(:)-test_true_y(:)).^2);

    test_frame_indices = find(is_test_frame);

    per_frame_errors(test_frame_indices) = euclidean_error;
    decoded_x_all(test_frame_indices) = decoded_x_position;
    decoded_y_all(test_frame_indices) = decoded_y_position;

end

per_frame_errors = single(per_frame_errors);

frames_with_valid_decoding = ~isnan(per_frame_errors);

median_error = median(per_frame_errors(frames_with_valid_decoding));
mean_error   = mean(per_frame_errors(frames_with_valid_decoding));
end


%% to concatanate as according to original method Hisa used
function combined_struct = safe_cat_structs(cell_array_of_structs)

has_content = ~cellfun(@isempty, cell_array_of_structs);
cell_array_of_structs = cell_array_of_structs(has_content);
number_of_entries = length(cell_array_of_structs);
if number_of_entries == 0, combined_struct = []; return; end
if number_of_entries == 1, combined_struct = cell_array_of_structs{1}; return; end

all_field_names = fieldnames(cell_array_of_structs{1});
for i = 2:number_of_entries
    all_field_names = union(all_field_names, fieldnames(cell_array_of_structs{i}), 'stable');
end
for i = 1:number_of_entries
    this_struct = cell_array_of_structs{i};
    for fi = 1:length(all_field_names)
        if ~isfield(this_struct, all_field_names{fi})
            this_struct.(all_field_names{fi}) = [];
        end
    end
    cell_array_of_structs{i} = orderfields(this_struct, all_field_names);
end
combined_struct = [cell_array_of_structs{:}];

end




%% Baysian using all cells (for 1st loop)
function [all_decoding_results, cached_session_data] = func_BayesDec(number_of_populations, number_of_sessions, results_directory, session_list, ...
    data_directory, minimum_speed_cm_per_s, arena_boundaries_cm, spatial_bin_size_cm, gaussian_kernel, number_of_cv_folds, ...
    minimum_cell_count, CellList_category_all, population_field_names, population_display_names, border_width_cm)
%%

all_decoding_results = cell(number_of_populations, 1);
for pop = 1:number_of_populations
    all_decoding_results{pop} = cell(number_of_sessions, 1);
end

% Cache preprocessed session data
cached_session_data = cell(number_of_sessions, 1);

scatter_plot_directory = fullfile(results_directory, 'ScatterPlots');
if ~exist(scatter_plot_directory, 'dir'), mkdir(scatter_plot_directory); end

total_timer = tic;

for session_idx = 1:number_of_sessions
    animal_name   = session_list(session_idx).animal_name;
    recording_day = session_list(session_idx).recording_day;
    mouse_index   = session_list(session_idx).mouse_index;
    group_label   = session_list(session_idx).group_label;
    this_filename = session_list(session_idx).filename;

    fprintf('--------------------------------------------------\n');
    fprintf('  Session %d/%d: %s day %d (%s)\n', ...
        session_idx, number_of_sessions, animal_name, recording_day, group_label);
    fprintf('--------------------------------------------------\n');

    try
        load_timer = tic;

        % Load calcium traces and tracking data
        loaded_data = load(fullfile(data_directory, this_filename));
        if ~isfield(loaded_data, 'CaTrace') || ~isfield(loaded_data, 'Table_video')
            fprintf('  ERROR: Missing CaTrace or Table_video\n'); continue;
        end

        calcium_timestamps  = single(loaded_data.CaTrace{:, 1});
        calcium_traces_raw  = single(loaded_data.CaTrace{:, 2:end});
        frame_interval_sec  = median(diff(calcium_timestamps));
        total_cells_in_file = size(calcium_traces_raw, 2);


        %Align tracking to calcium times;
        tracking_timestamps = loaded_data.Table_video.time_aligned;
        position_x_raw      = loaded_data.Table_video.x;
        position_y_raw      = loaded_data.Table_video.y;
        speed_1sec_mean     = loaded_data.Table_video.V_1secMean;
        [unique_tracking_times, unique_indices] = unique(tracking_timestamps, 'stable');
        position_x_aligned = interp1(unique_tracking_times, position_x_raw(unique_indices), calcium_timestamps, 'linear', NaN);
        position_y_aligned = interp1(unique_tracking_times, position_y_raw(unique_indices), calcium_timestamps, 'linear', NaN);
        speed_aligned      = interp1(unique_tracking_times, speed_1sec_mean(unique_indices), calcium_timestamps, 'linear', NaN);

        %  find thresholded motion frames
        is_moving_in_arena = speed_aligned >= minimum_speed_cm_per_s & ...
            ~isnan(position_x_aligned) & ~isnan(position_y_aligned) & ...
            position_x_aligned >= arena_boundaries_cm(1,1) & position_x_aligned <= arena_boundaries_cm(1,2) & ...
            position_y_aligned >= arena_boundaries_cm(2,1) & position_y_aligned <= arena_boundaries_cm(2,2);

        % make spatial bins
        x_bin_edges = arena_boundaries_cm(1,1):spatial_bin_size_cm:arena_boundaries_cm(1,2);
        y_bin_edges = arena_boundaries_cm(2,1):spatial_bin_size_cm:arena_boundaries_cm(2,2);
        number_of_x_bins = length(x_bin_edges) - 1;
        number_of_y_bins = length(y_bin_edges) - 1;
        x_bin_centers = x_bin_edges(1:end-1) + spatial_bin_size_cm / 2;
        y_bin_centers = y_bin_edges(1:end-1) + spatial_bin_size_cm / 2;
        [~, ~, x_bin_index_all_frames] = histcounts(position_x_aligned, x_bin_edges);
        [~, ~, y_bin_index_all_frames] = histcounts(position_y_aligned, y_bin_edges);

        %Using only valid movement frames with valid bin assignments
        movement_frame_indices = find(is_moving_in_arena);
        x_bin_during_movement  = x_bin_index_all_frames(movement_frame_indices);
        y_bin_during_movement  = y_bin_index_all_frames(movement_frame_indices);
        has_valid_bin = x_bin_during_movement > 0 & y_bin_during_movement > 0;
        movement_frame_indices = int32(movement_frame_indices(has_valid_bin));
        x_bin_during_movement  = x_bin_during_movement(has_valid_bin);
        y_bin_during_movement  = y_bin_during_movement(has_valid_bin);
        number_of_movement_frames = length(movement_frame_indices);

        true_position_x = x_bin_centers(x_bin_during_movement)';
        true_position_y = y_bin_centers(y_bin_during_movement)';

        % ==========================================================
        % Border vs Center classification
        % ==========================================================
        % border_width_cm = 15;
        arena_size_cm = arena_boundaries_cm(1,2);

        is_border_frame = ...
            true_position_x <= border_width_cm | ...
            true_position_x >= (arena_size_cm - border_width_cm) | ...
            true_position_y <= border_width_cm | ...
            true_position_y >= (arena_size_cm - border_width_cm);

        is_center_frame = ~is_border_frame;
        border_idx = find(is_border_frame); %index during movement
        center_idx = find(is_center_frame);

        % Safety check (avoid very small samples)
        % min_required_frames = 200;
        % if length(border_idx) < min_required_frames || length(center_idx) < min_required_frames
        %     warning('Too few border or center frames for reliable decoding.');
        % end

        %For cross-validation folds (contiguous temporal blocks, not random)
        frames_per_fold = floor(number_of_movement_frames / number_of_cv_folds);
        cv_fold_assignment = zeros(number_of_movement_frames, 1);
        for fold = 1:number_of_cv_folds
            if fold < number_of_cv_folds
                cv_fold_assignment((fold-1)*frames_per_fold+1 : fold*frames_per_fold) = fold;
            else
                cv_fold_assignment((fold-1)*frames_per_fold+1 : end) = fold;
            end
        end


        fprintf('Border frames: %d, Center frames: %d\n', length(border_idx), length(center_idx));
        fprintf('Folds in border: %s\n', mat2str(unique(cv_fold_assignment(border_idx))));
        fprintf('Folds in center: %s\n', mat2str(unique(cv_fold_assignment(center_idx))));

        %  Cache this session for subsampling
        this_session_cache = struct();
        this_session_cache.calcium_traces_raw    = calcium_traces_raw;
        this_session_cache.movement_frame_indices = movement_frame_indices;
        this_session_cache.x_bin_during_movement  = x_bin_during_movement;
        this_session_cache.y_bin_during_movement  = y_bin_during_movement;
        this_session_cache.true_position_x        = true_position_x;
        this_session_cache.true_position_y        = true_position_y;
        this_session_cache.cv_fold_assignment      = cv_fold_assignment;
        this_session_cache.frame_interval_sec      = frame_interval_sec;
        this_session_cache.number_of_x_bins        = number_of_x_bins;
        this_session_cache.number_of_y_bins        = number_of_y_bins;
        this_session_cache.x_bin_centers           = x_bin_centers;
        this_session_cache.y_bin_centers           = y_bin_centers;

        this_session_cache.is_border_frame           = is_border_frame; %index of border/center during running period

        cached_session_data{session_idx} = this_session_cache;

        fprintf('  Loaded in %.1f sec (%d movement frames)\n', toc(load_timer), number_of_movement_frames);

        % decode each cell population,
        cell_info_entry = CellList_category_all{mouse_index, recording_day};

        for pop = 1:number_of_populations
            % tic

            this_population_field = population_field_names{pop};
            this_population_label = population_display_names{pop};

            if ~isfield(cell_info_entry, this_population_field)
                fprintf('  [%s] Not present, skip\n', this_population_label); continue;
            end
            cell_table = cell_info_entry.(this_population_field);
            if height(cell_table) < minimum_cell_count
                fprintf('  [%s] Only %d cells, skip\n', this_population_label, height(cell_table)); continue;
            end

            selected_cell_ids = cell_table.CellInd;
            selected_cell_ids = selected_cell_ids(selected_cell_ids <= total_cells_in_file);
            number_of_selected_cells = length(selected_cell_ids);
            if number_of_selected_cells < minimum_cell_count
                fprintf('  [%s] %d valid cells < %d, skip\n', this_population_label, number_of_selected_cells, minimum_cell_count); continue;
            end

            % Z-score calcium traces 
            selected_traces = calcium_traces_raw(:, selected_cell_ids);
            trace_mean = mean(selected_traces, 1, 'omitnan');
            trace_std  = std(selected_traces, 0, 1, 'omitnan');
            trace_std(trace_std < 1e-10) = 1;
            zscored_traces = (selected_traces - trace_mean) ./ trace_std;
            zscored_traces_during_movement = zscored_traces(movement_frame_indices, :);


            % toc

            %%
            decoding_timer = tic;

            % ==========================================================
            % Decode using ALL frames (original behavior)
            % ==========================================================

            [median_decoding_error, mean_decoding_error, per_frame_errors, decoded_x, decoded_y, observed_x, observed_y] = ...
                run_bayesian_decoder(zscored_traces_during_movement, true_position_x, true_position_y, ...
                x_bin_during_movement, y_bin_during_movement, number_of_x_bins, number_of_y_bins, x_bin_centers, y_bin_centers, ...
                frame_interval_sec, gaussian_kernel, cv_fold_assignment, number_of_cv_folds);

            this_result = struct();
            this_result.animal_name          = animal_name;
            this_result.recording_day        = recording_day;
            this_result.group_label          = group_label;
            this_result.mouse_index          = mouse_index;
            this_result.population_label     = this_population_label;
            this_result.number_of_cells      = number_of_selected_cells;
            this_result.cell_ids             = selected_cell_ids;
            this_result.median_decoding_error = median_decoding_error;
            this_result.mean_decoding_error   = mean_decoding_error;
            this_result.per_frame_errors      = single(per_frame_errors);
            this_result.number_of_movement_frames = number_of_movement_frames;
            this_result.movement_frame_indices = movement_frame_indices;

            this_result.decoded_xy = single([decoded_x, decoded_y]);
            this_result.observed_xy = single([observed_x, observed_y]);
            this_result.unique_tracking_times = single(unique_tracking_times);
            this_result.unique_tracking_times_moving = single(unique_tracking_times(is_moving_in_arena));

            % % ==========================================================
            % % Decode using BORDER frames
            % % ==========================================================
            % 
            % [median_border, mean_border, err_border, dec_x_border, dec_y_border, obs_x_border, obs_y_border] = ...
            %     run_bayesian_decoder(zscored_traces_during_movement(border_idx,:), true_position_x(border_idx), true_position_y(border_idx), ...
            %     x_bin_during_movement(border_idx), y_bin_during_movement(border_idx), number_of_x_bins, number_of_y_bins, x_bin_centers, y_bin_centers, ...
            %     frame_interval_sec, gaussian_kernel, cv_fold_assignment(border_idx), number_of_cv_folds);
            % 
            % this_result.median_decoding_error_border = median_border;
            % this_result.mean_decoding_error_border   = mean_border;
            % this_result.per_frame_errors_border      = single(err_border);
            % 
            % % ==========================================================
            % % Decode using CENTER frames
            % % ==========================================================
            % 
            % [median_center, mean_center, err_center, dec_x_center, dec_y_center, obs_x_center, obs_y_center] = ...
            %     run_bayesian_decoder(zscored_traces_during_movement(center_idx,:), true_position_x(center_idx), ...
            %     true_position_y(center_idx), x_bin_during_movement(center_idx), y_bin_during_movement(center_idx), number_of_x_bins, number_of_y_bins, x_bin_centers, y_bin_centers, ...
            %     frame_interval_sec, gaussian_kernel, cv_fold_assignment(center_idx), number_of_cv_folds);
            % 
            % this_result.median_decoding_error_center = median_center;
            % this_result.mean_decoding_error_center   = mean_center;
            % this_result.per_frame_errors_center      = single(err_center);

            decoding_time = toc(decoding_timer);
            this_result.border_fraction = mean(is_border_frame);
            fprintf('  Decoding in %.1f sec\n', decoding_time);


            %% visualization

            % % Save observed-vs-decoded scatter plot ---(two plots,
            % % using scatter and density scatter
            % valid_decoded_frames = ~isnan(decoded_x) & ~isnan(observed_x);
            % if sum(valid_decoded_frames) > 10
            %     if strcmp(group_label, 'WildType')
            %         scatter_color = [0.2 0.4 0.8];
            %     else
            %         scatter_color = [0.8 0.2 0.2];
            %     end
            %
            %
            %     density_figure = figure('Visible', 'off', 'Position', [100 100 900 400]);
            %     number_of_density_bins = 40;
            %     density_edges = linspace(0, 100, number_of_density_bins + 1);
            %
            %     subplot(1,2,1);
            %     density_counts_x = histcounts2(...
            %         observed_x(valid_decoded_frames), decoded_x(valid_decoded_frames), ...
            %         density_edges, density_edges);
            %     imagesc(density_edges(1:end-1) + 1.25, density_edges(1:end-1) + 1.25, density_counts_x');
            %     axis xy square; xlim([0 100]); ylim([0 100]);
            %     hold on; plot([0 100], [0 100], 'w--', 'LineWidth', 1.5);
            %     xlabel('Observed X (cm)'); ylabel('Decoded X (cm)');
            %     colormap(gca, hot); colorbar; caxis([0 prctile(density_counts_x(:), 99)]);
            %     title('X density');
            %
            %     subplot(1,2,2);
            %     density_counts_y = histcounts2(...
            %         observed_y(valid_decoded_frames), decoded_y(valid_decoded_frames), ...
            %         density_edges, density_edges);
            %     imagesc(density_edges(1:end-1) + 1.25, density_edges(1:end-1) + 1.25, density_counts_y');
            %     axis xy square; xlim([0 100]); ylim([0 100]);
            %     hold on; plot([0 100], [0 100], 'w--', 'LineWidth', 1.5);
            %     xlabel('Observed Y (cm)'); ylabel('Decoded Y (cm)');
            %     colormap(gca, hot); colorbar; caxis([0 prctile(density_counts_y(:), 99)]);
            %     title('Y density');
            %
            %     sgtitle(sprintf('%s | %s day %d | %s | n=%d cells | error=%.1f cm', ...
            %         this_population_label, animal_name, recording_day, group_label, ...
            %         number_of_selected_cells, median_decoding_error), 'FontSize', 11);
            %     saveas(density_figure, fullfile(scatter_plot_directory, ...
            %         sprintf('Density_%s_%s_d%d.png', this_population_label, animal_name, recording_day)));
            %     close(density_figure);
            % end
            %


            all_decoding_results{pop}{session_idx} = this_result;
        end

    catch error_info
        fprintf('  *** ERROR: %s\n', error_info.message);
        fprintf('  *** at %s line %d\n', error_info.stack(1).name, error_info.stack(1).line);
    end
end
fprintf('\nAll sessions decoded in %.0f sec\n', toc(total_timer));

end
