function [HD_Score, moveHD_Score]=HD_score_nofilter_RL(DIR,CDir,SampleName, t)

cd(DIR);

load(fullfile(DIR,'Angle_FR_Score_nofilter.mat'), ...
    "Bined_HeadBodyDirection_FR", ...
    "HD_Score_Bin6", ...
    "moveBined_HeadBodyDirection_FR", ...
    "moveHD_Score_Bin6");

load(fullfile(DIR,'ST_dF_grid_aut_data.mat'), ...
    "Original_Cell_ID");

numCells = size(HD_Score_Bin6,1);
HD_Score = nan(numCells,6);
moveHD_Score = nan(numCells,6);

ang_sd_deg = nan(numCells,1);
move_ang_sd_deg = nan(numCells,1);
body_ang_sd_deg = nan(numCells,1);
body_move_ang_sd_deg = nan(numCells,1);
HDcell_num = 1;

for c = 1:numCells
    % Columns: InfoRate, Watson's U2-test p-value, Watson's U2 score,
    % Rayleigh p-value, Rayleigh z-score, weighted angular SD.
    HD_Score(c,1:5) = HD_Score_Bin6{c,1};
    moveHD_Score(c,1:5) = moveHD_Score_Bin6{c,1};

    direction_firing_matrix = Bined_HeadBodyDirection_FR{c,1}(:,1:2);
    ang_sd_deg(c) = weighted_angular_std_deg(direction_firing_matrix);
    HD_Score(c,6) = ang_sd_deg(c);

    move_direction_firing_matrix = moveBined_HeadBodyDirection_FR{c,1}(:,1:2);
    move_ang_sd_deg(c) = weighted_angular_std_deg(move_direction_firing_matrix);
    moveHD_Score(c,6) = move_ang_sd_deg(c);

    direction_firing_matrix = Bined_HeadBodyDirection_FR{c,2}(:,1:2);
    body_ang_sd_deg(c) = weighted_angular_std_deg(direction_firing_matrix);

    move_direction_firing_matrix = moveBined_HeadBodyDirection_FR{c,2}(:,1:2);
    body_move_ang_sd_deg(c) = weighted_angular_std_deg(move_direction_firing_matrix);

    if (HD_Score(c,2) < 0.05 && HD_Score(c,4) < 0.05 && HD_Score(c,6) < 40)
        oriC = Original_Cell_ID(c);
        clf
        figure('Position', [0, 0, 800, 900]);

        HeadTitle = "Head direction FR (Bin6)";
        ax = subplot(4,4,[1 2 5 6], polaraxes);
        plot_polar_spike4(Bined_HeadBodyDirection_FR{c,1}(:,[1 2]),oriC,HD_Score_Bin6{c,1},ang_sd_deg(c),HeadTitle,ax);

        HeadTitle = "Head direction FR (Bin6) >2cm/s";
        ax = subplot(4,4,[3 4 7 8], polaraxes);
        plot_polar_spike4(moveBined_HeadBodyDirection_FR{c,1}(:,[1 2]),oriC,moveHD_Score_Bin6{c,1},move_ang_sd_deg(c),HeadTitle,ax);

        HeadTitle = "Body direction FR (Bin6)";
        ax = subplot(4,4,[9 10 13 14], polaraxes);
        plot_polar_spike4(Bined_HeadBodyDirection_FR{c,2}(:,[1 2]),oriC,HD_Score_Bin6{c,2},body_ang_sd_deg(c),HeadTitle,ax);

        HeadTitle = "Body direction FR (Bin6) >2cm/s";
        ax = subplot(4,4,[11 12 15 16], polaraxes);
        plot_polar_spike4(moveBined_HeadBodyDirection_FR{c,2}(:,[1 2]),oriC,moveHD_Score_Bin6{c,2},body_move_ang_sd_deg(c),HeadTitle,ax);

        nameTitle = strcat(CDir,"\nofilter_HDcell\",SampleName, "_T",num2str(t)," cell#", num2str(c), ".jpg");
        exportgraphics(gcf,nameTitle)

        HDcell_num = HDcell_num + 1;
    end
end

HDcell_num
numCells
end
