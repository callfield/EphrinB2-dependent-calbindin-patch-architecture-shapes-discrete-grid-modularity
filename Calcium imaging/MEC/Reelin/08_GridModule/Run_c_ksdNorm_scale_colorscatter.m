clear all; close all;

addpath('function')
load("../Data.mat")
load("GridMod.mat")
CDir = pwd;

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

%% Gmod trial grid with KSD regression slope
gp1_Mod{5,3} = [];
gp1 = [];
gp1_1 = [];
gp1_2 = [];
for s = 1:5
    for t = 1:3
        if isempty(gp1_Mod{s,t}) == 0
            data = GROUP1{s,t}(GROUP1{s,t}(:,3) > 0, 1:13);

            Mod1 = gp1_Mod{s,t}{1};
            Mod2 = gp1_Mod{s,t}{2};
            peaks_GM1scale = gp1_GSpeak{s,t}(1,1); % Peak of Grid scale of Gmod1.
            norm_data = data;
            norm_data(:,3) = data(:,3) ./ peaks_GM1scale;

            gp1 = [gp1; norm_data];
            gp1_1 = [gp1_1; norm_data(Mod1,:)];
            gp1_2 = [gp1_2; norm_data(Mod2,:)];
        end
    end
end

gp2 = [];
gp2_1 = [];
gp2_2 = [];
for s = 1:7
    for t = 1:3
        if isempty(gp2_Mod{s,t}) == 0
            data = GROUP2{s,t}(GROUP2{s,t}(:,3) > 0, 1:13);

            Mod1 = gp2_Mod{s,t}{1};
            Mod2 = gp2_Mod{s,t}{2};

            peaks_GM1scale = gp2_GSpeak{s,t}(1,1); % Peak of Grid scale of Gmod1.
            norm_data = data;
            norm_data(:,3) = data(:,3) ./ peaks_GM1scale;

            gp2 = [gp2; norm_data];
            gp2_1 = [gp2_1; norm_data(Mod1,:)];
            gp2_2 = [gp2_2; norm_data(Mod2,:)];
        end
    end
end

mkdir NormGS_Gmodtrial

% Data.
data1 = gp1(:,3);
data2 = gp2(:,3);

% Set a common x-axis range for both groups.
xmin = min([data1; data2]);
xmax = max([data1; data2]);
x_values = linspace(xmin, xmax, 200);

% Compute KDE.
[f1, xi1] = ksdensity(data1, x_values, 'Bandwidth', 0.05);
[f2, xi2] = ksdensity(data2, x_values, 'Bandwidth', 0.05);

% Plot.
figure('Position', [100, 100, 600, 300]);
hold on;
area(xi1, f1, 'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 0.5);
plot(xi1, f1, 'Color', [0.3 0.3 0.3] * 0.8, 'LineWidth', 2);
area(xi2, f2, 'FaceColor', [0.9 0.3 0.3], 'FaceAlpha', 0.5);
plot(xi2, f2, 'Color', [0.9 0.3 0.3] * 0.8, 'LineWidth', 2);

xlabel('Grid Scale');
ylabel('Density');
[h, p, ksstat] = kstest2(data1, data2);

% Show the p-value and test result in the title.
title(sprintf('Kolmogorov-Smirnov test: p = %.4f, h = %d', p, h));

exportgraphics(gcf, strcat("NormGS_Gmodtrial\NormGS.pdf"))
close all

xlims = [0.5 2.5];
ylims = [0 800];
peaks = ksdnorm200(gp1);
plotColor_GridScaleVsLocation_Mod1norm(gp1, peaks, gp1_1, gp1_2, xlims, ylims);
exportgraphics(gcf, strcat("NormGS_Gmodtrial\NormGS ScatterwithSlope group1.pdf"))

peaks = ksdnorm200(gp2);
plotColor_GridScaleVsLocation_Mod1norm(gp2, peaks, gp2_1, gp2_2, xlims, ylims);
exportgraphics(gcf, strcat("NormGS_Gmodtrial\NormGS ScatterwithSlope group2.pdf"))

close all
