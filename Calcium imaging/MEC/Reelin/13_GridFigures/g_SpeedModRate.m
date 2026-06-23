close all; clear all


load('../08_GridModule/GridMod.mat');
load("../Data.mat")
load("../09_Grid_module_speed/ZSpeedHz_each1.mat") 

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

%%

gp1_ZSpeedHz_amean=cell(3);
gp1_cell_profiles=cell(3);
ii=1;
for s=1:5 %Animal
    for t=1:3 %Trials
        gcell= find(GROUP1{s,t}(:,4)~=0);
        for m=1

            tmp=gp1_ZSpeedHz{s,t}(1:20, gcell);
            tmp(tmp==0)=NaN;% convert zero to nan before mean
            tmp2(1,1)=mean(tmp(1:5,:),'all','omitmissing');
            tmp2(2,1)=mean(tmp(6:10,:),'all','omitmissing');
            tmp2(3,1)=mean(tmp(11:15,:),'all','omitmissing');
            tmp2(4,1)=mean(tmp(16:20,:),'all','omitmissing');

            gp1_ZSpeedHz_amean{m}=[gp1_ZSpeedHz_amean{m}, tmp2];
            tmp2=[];

            gp1_cell_profiles{m} = [gp1_cell_profiles{m}, tmp];

        end
    end
end

gp2_ZSpeedHz_amean=cell(2,1);
gp2_cell_profiles = cell(2,1);
ii=1;
for s=1:7
    for t=1:3
        if isempty(GROUP2{s,t}) == 0
            gcell= find(GROUP2{s,t}(:,4)~=0);
            for m=1 %1:2

                tmp=gp2_ZSpeedHz{s,t}(1:20, gcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean
                tmp2(1,1)=mean(tmp(1:5,:),'all','omitmissing');
                tmp2(2,1)=mean(tmp(6:10,:),'all','omitmissing');
                tmp2(3,1)=mean(tmp(11:15,:),'all','omitmissing');
                tmp2(4,1)=mean(tmp(16:20,:),'all','omitmissing');
                gp2_ZSpeedHz_amean{m}=[gp2_ZSpeedHz_amean{m}, tmp2];

                tmp2=[];
                gp2_cell_profiles{m} = [gp2_cell_profiles{m}, tmp];
            end

        end
    end
end


%%

close all

data1_gp1 = gp1_ZSpeedHz_amean{1,1};
data2_gp1 = gp1_ZSpeedHz_amean{2,1};


%% percentage of speed modulated cell
% =========================================================
% Session-based analysis
% each trial is treated as independent
% ==========================================================

close all

alpha = 0.05;

% ---------- Group1 sessions ----------
gp1_session_frac = [];
for s = 1:5
    for t = 1:3
        if isempty(gp1_ZSpeedHz{s,t})
            continue
        end
        tmp = gp1_ZSpeedHz{s,t}(1:20,:);
        [corrs, pvals] = compute_speed_corr_single(tmp);
        is_mod = pvals < alpha;
        frac_mod = sum(is_mod) / length(is_mod);
        gp1_session_frac(end+1) = frac_mod;
    end
end


% ---------- Group2 sessions ----------
gp2_session_frac = [];
for s = 1:7
    for t = 1:3
        if isempty(gp2_ZSpeedHz{s,t})
            continue
        end
        tmp = gp2_ZSpeedHz{s,t}(1:20,:);
        [corrs, pvals] = compute_speed_corr_single(tmp);
        is_mod = pvals < alpha;
        frac_mod = sum(is_mod) / length(is_mod);
        gp2_session_frac(end+1) = frac_mod;
    end
end


% ---------- statistics ----------

[p_mod,~,stats] = ranksum( ...
    gp1_session_frac, ...
    gp2_session_frac);

fprintf('\n');
fprintf('Session-based comparison\n');

fprintf('Group1 = %.3f ± %.3f\n', ...
    mean(gp1_session_frac), ...
    std(gp1_session_frac)/sqrt(length(gp1_session_frac)));

fprintf('Group2    = %.3f ± %.3f\n', ...
    mean(gp2_session_frac), ...
    std(gp2_session_frac)/sqrt(length(gp2_session_frac)));

fprintf('Ranksum p = %.4f\n', p_mod);


% =========================================================
% Figure
% ==========================================================

close all

CGroup1 = [1 1 1];
CGroup2 = [255 168 60]/255;

fig = figure('Position',[300 400 150 180]);

hold on

x = [0.9 2.1];

means = [
    mean(gp1_session_frac)
    mean(gp2_session_frac)
    ];

b = bar(x, means, ...
    'FaceColor','flat', ...
    'BarWidth',0.6);

b.CData(1,:) = CGroup1;
b.CData(2,:) = CGroup2;


% error bar

sem1 = std(gp1_session_frac) ...
    / sqrt(length(gp1_session_frac));

sem2 = std(gp2_session_frac) ...
    / sqrt(length(gp2_session_frac));

errorbar(x, means, ...
    [sem1; sem2], ...
    'k', ...
    'LineStyle','none', ...
    'LineWidth',1, 'CapSize', 15);


% scatter
ytol = 0.01;
dx = 0.10;
x1 = SimpleBeeSwarm(gp1_session_frac, x(1), ytol, dx);
x2 = SimpleBeeSwarm(gp2_session_frac, x(2), ytol, dx);
scatter(x1, gp1_session_frac, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
scatter(x2, gp2_session_frac, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)



% p value
y = max([gp1_session_frac gp2_session_frac]) + 0.08;
% plot([x(1) x(2)], [y y], 'k')
% text(mean(x), y+0.03, ...
%     ['p = ' num2str(p_mod,'%.3f')], ...
%     'HorizontalAlignment','center')



% ylabel('Fraction')
xticks(x)
% xticklabels({'Group1','Group2'})

ylim([0 1])
xlim([0 3])

yticks([0 0.5 1])

xticklabels([])
yticklabels([])

box off
set(gca,'LineWidth',1)

s = strcat("SpeedModRatio_Group2.emf");
exportgraphics(fig,  s);



A = gp1_session_frac;
B = gp2_session_frac;
mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
fprintf('Group2: %.3f ± %.3f\n', mB, sdB);


%% ========== statistics ==========

% NaN除去
Y1 = gp1_session_frac(~isnan(gp1_session_frac));
Y2 = gp2_session_frac(~isnan(gp2_session_frac));

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

% smaller U
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








%%
function [corrs, pvals] = compute_speed_corr_single(tmp)

speed_axis = (1:20)';

nCells = size(tmp,2);

corrs = nan(1,nCells);
pvals = nan(1,nCells);

for i = 1:nCells

    y = tmp(:,i);

    valid = ~isnan(y);

    if sum(valid) < 5
        continue
    end

    [r,p] = corr(speed_axis(valid), y(valid), ...
        'type','Spearman');

    corrs(i) = r;
    pvals(i) = p;

end

end



%%
function xj = SimpleBeeSwarm(y, xcenter, ytol, dx)

% y       : y values
% xcenter : bar center (e.g. 1 or 2)
% ytol    : "same height" threshold
% dx      : horizontal spacing

y = y(:);

xj = zeros(size(y));

% sort by y
[ys, idx] = sort(y);

groups = {};
g = 1;
groups{g} = idx(1);

for i = 2:length(ys)

    if abs(ys(i) - ys(i-1)) < ytol
        groups{g}(end+1) = idx(i);
    else
        g = g + 1;
        groups{g} = idx(i);
    end
end

% assign x offsets
for g = 1:length(groups)

    ids = groups{g};
    n = length(ids);

    if n == 1
        offsets = 0;
    else
        offsets = ((1:n) - mean(1:n)) * dx;
    end

    xj(ids) = xcenter + offsets;
end

end