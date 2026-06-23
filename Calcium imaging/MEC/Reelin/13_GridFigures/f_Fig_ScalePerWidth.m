%% ========== Load data ==========
clear; close all hidden


pos = [300 500 130 160];


%% Group2

Ccont = [1 1 1];
CGroup2 = [255 168 60]/255;
% Ccont = [1 1 1];
% CGroup2 = [240 134 134]/255;


% Gscale_width_ratio_all.csv: 1  Gscale/gwidth of all grid cells in Group1, 2  Gscale/gwidth of all grid cells in Group2 
Data_temp = readmatrix('Gscale_width_ratio_all.csv');

X = Data_temp(:,1);
Y = Data_temp(:,2);

ydata = [X; Y];

group = [ ...
    repmat("Group1", length(X), 1); ...
    repmat("Group2", length(Y), 1)];

group = categorical(group);


valid = ~isnan(ydata);
ydata = ydata(valid);
group = group(valid);


fig2 = figure('Position', pos);

Y = {X, Y};

hold on

% vwidth = 0.3;

xgroupdata = categorical(repelem(["group1";"group2"],[length(Y{1}), length(Y{2})]));
% v = violinplot_2(cell2mat(Y'), xgroupdata, 'ShowData', false, 'ViolinColor', [1 1 1]*0.8,...
%     'BoxColor', [1 1 1]*0.4, 'Width', 0.3, 'EdgeColor', [1 1 1]*1, 'ViolinAlpha', 1);

v = violinplot_2(cell2mat(Y'), xgroupdata, 'ShowData', false, 'ViolinColor', [1 1 1]*0.8,...
    'BoxColor', [1 1 1]*0.4, 'Width', 0.3, 'EdgeColor', [1 1 1]*1, 'ViolinAlpha', 1, 'BandWidth', 0.15);

v(2).ViolinColor = {CGroup2};

ax = gca;
xlim([0.4 2.6])
box off
xticks('')

ylim([0 4])
yticks([0 2 4])
yticklabels({})


% ========== statistics ==========
% ranksum
[p_rank, h_rank, stats_rank] = ranksum(Y{1}, Y{2});

% effect size for ranksum (rank-biserial style z/sqrt(N))
effectsize_rank = stats_rank.zval / sqrt(length(group));

% t-test
[~, p_ttest, ci_ttest, stats_ttest] = ttest2(Y{1}, Y{2});

% Cohen's d
n1 = length(Y{1});
n2 = length(Y{2});

mean1 = mean(Y{1}, 'omitnan');
mean2 = mean(Y{2}, 'omitnan');

sd1 = std(Y{1}, 'omitnan');
sd2 = std(Y{2}, 'omitnan');

% pooled SD
spooled = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2));
cohen_d = (mean1 - mean2) / spooled;

% ========== print results ==========
fprintf('===== t-test =====\n');
fprintf('t(%d) = %.4f\n', stats_ttest.df, stats_ttest.tstat);
fprintf('p = %.6g\n', p_ttest);
fprintf('95%% CI = [%.4f, %.4f]\n\n', ci_ttest(1), ci_ttest(2));

fprintf('===== ranksum =====\n');
fprintf('z = %.4f\n', stats_rank.zval);
fprintf('p = %.6g\n', p_rank);
fprintf('effect size (z/sqrt(N)) = %.4f\n\n', effectsize_rank);

fprintf('===== Cohen''s d =====\n');
fprintf('Cohen''s d = %.4f\n', cohen_d);

s = strcat("FigS GridScalePerWidth_Group2.emf");
exportgraphics(fig2, s);


A = Y{1};
    B = Y{2};
    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    sdA = std(A,'omitmissing');
    sdB = std(B,'omitmissing');
    fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
    fprintf('Group2: %.3f ± %.3f\n', mB, sdB);


% ========== statistics2 ==========
Y1 = Y{1}(~isnan(Y{1}));
Y2 = Y{2}(~isnan(Y{2}));

% sample size
n1 = length(Y1);
n2 = length(Y2);

fprintf('===== Sample size =====\n');
fprintf('Group1 n = %d\n', n1);
fprintf('Group2    n = %d\n\n', n2);

% =========================================================
% Mann-Whitney U test (ranksum)
% =========================================================

[p_rank, h_rank, stats_rank] = ranksum(Y1, Y2);

% MATLAB ranksum returns rank sum
% convert to Mann-Whitney U
R1 = stats_rank.ranksum;
U1 = R1 - n1*(n1+1)/2;
U2 = n1*n2 - U1;
U = min(U1, U2);

% effect size
effectsize_rank = stats_rank.zval / sqrt(n1+n2);

fprintf('===== Mann-Whitney U test =====\n');
fprintf('U = %.4f\n', U);
fprintf('z = %.4f\n', stats_rank.zval);
fprintf('p = %.6g\n', p_rank);
fprintf('effect size (z/sqrt(N)) = %.4f\n\n', effectsize_rank);

% =========================================================
% t-test
% =========================================================
[~, p_ttest, ci_ttest, stats_ttest] = ttest2(Y1, Y2);
fprintf('===== t-test =====\n');
fprintf('t(%d) = %.4f\n', round(stats_ttest.df), stats_ttest.tstat);
fprintf('p = %.6g\n', p_ttest);
fprintf('95%% CI = [%.4f, %.4f]\n\n', ci_ttest(1), ci_ttest(2));

% =========================================================
% Cohen's d
% =========================================================
mean1 = mean(Y1, 'omitnan');
mean2 = mean(Y2, 'omitnan');
sd1 = std(Y1, 'omitnan');
sd2 = std(Y2, 'omitnan');

% pooled SD
spooled = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2));
cohen_d = (mean1 - mean2) / spooled;

fprintf('===== Cohen''s d =====\n');
fprintf('Cohen''s d = %.4f\n', cohen_d);



