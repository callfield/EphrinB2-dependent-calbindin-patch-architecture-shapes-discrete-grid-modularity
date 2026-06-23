clearvars;

scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'function'))

load(fullfile(scriptDir, "..", "Data.mat"))

% GROUP1, GROUP2:
% 1: DV position (um)
% 2: ML position (um)
% 3: Grid Scale
% 4: Grid Orientation
% 5: Grid Width
% 6:8 Three Orientation values
% 9:11 Three Grid Scale values
% 12: Grid Score
% 13: z Grid Score

%% Grid orientation scatter
% Keep grid cells with positive grid scale, then analyze column 4.
gp1Parts = cell(7, 3);
for s = 1:7
    for t = 1:3
        trialData = GROUP1{s,t};
        if isempty(trialData)
            gp1Parts{s,t} = zeros(0, 13);
        else
            gp1Parts{s,t} = trialData(trialData(:,3) > 0, 1:13);
        end
    end
end
gp1 = vertcat(gp1Parts{:});

gp2Parts = cell(7, 3);
for s = 1:7
    for t = 1:3
        trialData = GROUP2{s,t};
        if isempty(trialData)
            gp2Parts{s,t} = zeros(0, 13);
        else
            gp2Parts{s,t} = trialData(trialData(:,3) > 0, 1:13);
        end
    end
end
gp2 = vertcat(gp2Parts{:});

data1 = gp1;
data2 = gp2;

close all
figure('Position', [100, 100, 250, 500]);
xlims = [-30 30];
ylims = [0 800];

plot_GridOriVsLocation(data1, xlims, ylims);
exportgraphics(gcf, fullfile(scriptDir, "G_ori_scatter_group1_all_trial.pdf"))
close all

figure('Position', [100, 100, 250, 500]);
plot_GridOriVsLocation(data2, xlims, ylims);
exportgraphics(gcf, fullfile(scriptDir, "G_ori_scatter_group2_all_trial.pdf"))
close all

%% Kernel-smoothed grid orientation distribution
bandwid = 3.5;
numSamples = 1000;
rng(1);  % Reproducible random samples from the fitted kernel distributions.

ori1 = data1(:,4);
ori2 = data2(:,4);

ori1 = ori1(isfinite(ori1));
ori2 = ori2(isfinite(ori2));

pd1 = fitdist(ori1, 'Kernel', 'BandWidth', bandwid);
normOri1 = [ori1; random(pd1, numSamples, 1)];

pd2 = fitdist(ori2, 'Kernel', 'BandWidth', bandwid);
normOri2 = [ori2; random(pd2, numSamples, 1)];

clf
figure('Position', [100, 100, 500, 200]);
hold on
plotDistribution_ksd(normOri1, bandwid, [0.3 0.3 0.3]);
plotDistribution_ksd(normOri2, bandwid, [0.9 0.3 0.3]);

xlim([-30 30])
ax = gca;
ax.XAxis.FontSize = 11;
yticks([])

title("All trial", 'FontSize', 14);
xlabel("Grid Orientation (deg)", 'FontSize', 18);
ylabel("KSD probability", 'FontSize', 18);

exportgraphics(gcf, fullfile(scriptDir, "G_ori_KSD_all_trial.pdf"));
clf

[H, p, ksstat] = kstest2(normOri1, normOri2);
fprintf('KSD-sampled KS test: H = %d, p = %.4f, KS stat = %.4f\n', H, p, ksstat)
