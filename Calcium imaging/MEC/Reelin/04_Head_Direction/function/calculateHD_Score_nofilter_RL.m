function [HeadBodyDirection_FR, Bined_HeadBodyDirection_FR, HD_Score_Bin1, HD_Score_Bin6] = calculateHD_Score_nofilter_RL(DIR, caFr)
    % Move to the specified directory.
    cd(DIR);
    mkdir("Direction/Head_direction_nofilter");
    % mkdir("Direction/Body_direction");

    % Load angle data.
    load("Angle.mat", "HeadAngle", "BodyAngle");

    % Load recording data and metadata.
    load("ST_dF_grid_aut_data.mat", "lk", "Original_Cell_ID");
    Trk = csvread("ST_PCI_Ca_behav_track.csv");   % Tracking data such as position and speed
    dF  = csvread("ST_PCI_noDup_dF.csv");         % De-duplicated calcium traces
    mFrame = find(Trk(:,4) >= 2);                  % Frames with speed >= 2 cm/s

    [frameNum, cellNum] = size(dF);
    HeadBodyDirection_FR       = cell(cellNum, 2);
    Bined_HeadBodyDirection_FR = cell(cellNum, 2);
    HD_Score_Bin1              = cell(cellNum, 2);
    HD_Score_Bin6              = cell(cellNum, 2);

    moveHeadBodyDirection_FR       = cell(cellNum, 2);
    moveBined_HeadBodyDirection_FR = cell(cellNum, 2);
    moveHD_Score_Bin1              = cell(cellNum, 2);
    moveHD_Score_Bin6              = cell(cellNum, 2);

    ang_sd_deg = nan(cellNum, 1);
    move_ang_sd_deg = nan(cellNum, 1);
    body_ang_sd_deg = nan(cellNum, 1);
    body_move_ang_sd_deg = nan(cellNum, 1);

    for c = 1:cellNum
        % Build a binary spike train at the imaging frame rate.
        SpikeFrame = round(lk{c} * caFr);
        tmpSpike   = zeros(frameNum, 1);
        tmpSpike(SpikeFrame) = 1;

        % Optional spike smoothing is disabled here to keep raw spike events.
        % before = ceil(caFr/3); after = ceil(caFr/3);
        % filtered_tmpSpike = gaussian_filter_custom(tmpSpike, before, after);
        filtered_tmpSpike = tmpSpike;

        % Compute per-angle firing and HD scores (Bin1 and Bin6): head direction.
        DATA = [round(HeadAngle), filtered_tmpSpike];
        [HeadBodyDirection_FR{c,1}, Bined_HeadBodyDirection_FR{c,1}, HD_Score_Bin1{c,1}, HD_Score_Bin6{c,1}] = ...
            compute_angle_spike_average_RL(DATA, caFr);

        % Compute per-angle firing and HD scores: body direction.
        DATA = [round(BodyAngle), filtered_tmpSpike];
        [HeadBodyDirection_FR{c,2}, Bined_HeadBodyDirection_FR{c,2}, HD_Score_Bin1{c,2}, HD_Score_Bin6{c,2}] = ...
            compute_angle_spike_average_RL(DATA, caFr);

        % Only moving frames (speed >= 2 cm/s): head direction.
        DATA = [round(HeadAngle(mFrame,:)), filtered_tmpSpike(mFrame)];
        [moveHeadBodyDirection_FR{c,1}, moveBined_HeadBodyDirection_FR{c,1}, moveHD_Score_Bin1{c,1}, moveHD_Score_Bin6{c,1}] = ...
            compute_angle_spike_average_RL(DATA, caFr);

        % Only moving frames: body direction.
        DATA = [round(BodyAngle(mFrame,:)), filtered_tmpSpike(mFrame)];
        [moveHeadBodyDirection_FR{c,2}, moveBined_HeadBodyDirection_FR{c,2}, moveHD_Score_Bin1{c,2}, moveHD_Score_Bin6{c,2}] = ...
            compute_angle_spike_average_RL(DATA, caFr);

        % Weighted angular SD from binned head-direction firing rates.
        direction_firing_matrix = Bined_HeadBodyDirection_FR{c,1}(:, 1:2);
        ang_sd_deg(c) = weighted_angular_std_deg(direction_firing_matrix);

        % Weighted angular SD for moving frames only.
        move_direction_firing_matrix = moveBined_HeadBodyDirection_FR{c,1}(:, 1:2);
        move_ang_sd_deg(c) = weighted_angular_std_deg(move_direction_firing_matrix);

        % Weighted angular SD for body direction.
        direction_firing_matrix = Bined_HeadBodyDirection_FR{c,2}(:, 1:2);
        body_ang_sd_deg(c) = weighted_angular_std_deg(direction_firing_matrix);

        % Weighted angular SD for body direction during moving frames.
        move_direction_firing_matrix = moveBined_HeadBodyDirection_FR{c,2}(:, 1:2);
        body_move_ang_sd_deg(c) = weighted_angular_std_deg(move_direction_firing_matrix);

        % Visualization.
        oriC = Original_Cell_ID(c);
        clf;
        figure('Position', [0, 0, 800, 900]);  % [left, bottom, width, height]

        HeadTitle = "Head direction FR (Bin6)";
        ax = subplot(4, 4, [1 2 5 6], polaraxes);
        plot_polar_spike4(Bined_HeadBodyDirection_FR{c,1}(:, [1 2]), oriC, HD_Score_Bin6{c,1}, ang_sd_deg(c), HeadTitle, ax);

        HeadTitle = "Head direction FR (Bin6) >2 cm/s";
        ax = subplot(4, 4, [3 4 7 8], polaraxes);
        plot_polar_spike4(moveBined_HeadBodyDirection_FR{c,1}(:, [1 2]), oriC, moveHD_Score_Bin6{c,1}, move_ang_sd_deg(c), HeadTitle, ax);

        HeadTitle = "Body direction FR (Bin6)";
        ax = subplot(4, 4, [9 10 13 14], polaraxes);
        plot_polar_spike4(Bined_HeadBodyDirection_FR{c,2}(:, [1 2]), oriC, HD_Score_Bin6{c,2}, body_ang_sd_deg(c), HeadTitle, ax);

        HeadTitle = "Body direction FR (Bin6) >2 cm/s";
        ax = subplot(4, 4, [11 12 15 16], polaraxes);
        plot_polar_spike4(moveBined_HeadBodyDirection_FR{c,2}(:, [1 2]), oriC, moveHD_Score_Bin6{c,2}, body_move_ang_sd_deg(c), HeadTitle, ax);

        nameTitle = strcat("Direction/Head_direction_nofilter/Direction_score cell#", num2str(c), ".jpg");
        exportgraphics(gcf, nameTitle);

        % Periodically close figures to free memory.
        if mod(c, 30) == 0
            close all;
        end
    end

    close all;

    % Save full outputs, including no-speed-filter and moving-only variants.
    save("Angle_FR_Score_nofilter.mat", "HeadBodyDirection_FR", "Bined_HeadBodyDirection_FR", ...
        "HD_Score_Bin1", "HD_Score_Bin6", ...
        "moveHeadBodyDirection_FR", "moveBined_HeadBodyDirection_FR", ...
        "moveHD_Score_Bin1", "moveHD_Score_Bin6");
end
