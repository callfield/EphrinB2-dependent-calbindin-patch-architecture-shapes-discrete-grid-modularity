%% ========== Load data ==========
clear; close all hidden

GROUP1 = load("../Group1_CellTables.mat");
GROUP2 = load("../Group2_CellTables.mat");
GROUP1_T = GROUP1.GROUP1_T;
GROUP2_T = GROUP2.GROUP2_T;


min_modulenum = 2;
mode = 'adjacent'; % 'all' : use only module 1 and 2.


%% Use only module 1 and 2


for m = 1:numel(GROUP1_T)
    if ~isempty(GROUP1_T{m})
        idx3 = GROUP1_T{m}.Grid_module == 3;
        GROUP1_T{m}.Grid_module(idx3) = 0;
    end
end

for m = 1:numel(GROUP2_T)
    if ~isempty(GROUP2_T{m})
        idx3 = GROUP2_T{m}.Grid_module == 3;
        GROUP2_T{m}.Grid_module(idx3) = 0;
    end
end

GROUP1_T = GROUP1_T;
GROUP2_T = GROUP2_T;


%% max grid module number
MaxModule_Group1 = nan(5,3);
Num_Group1_M1 = nan(5,3);
Num_Group1_M2 = nan(5,3);
Num_Group1_M3 = nan(5,3);
for MouseI = 1:5
    for Trial_I = 1:3
        MaxModule_Group1(MouseI,Trial_I) = max(GROUP1_T{MouseI,Trial_I}.Grid_module);
        Num_Group1_M1(MouseI,Trial_I) = length(find(GROUP1_T{MouseI,Trial_I}.Grid_module == 1));
        Num_Group1_M2(MouseI,Trial_I) = length(find(GROUP1_T{MouseI,Trial_I}.Grid_module == 2));
        Num_Group1_M3(MouseI,Trial_I) = length(find(GROUP1_T{MouseI,Trial_I}.Grid_module == 3));
    end
end

MaxModule_Gp2 = nan(7,3);
Num_Gp2_M1 = nan(7,3);
Num_Gp2_M2 = nan(7,3);
Num_Gp2_M3 = nan(7,3);
for MouseI = 1:7
    for Trial_I = 1:3
        try
            MaxModule_Gp2(MouseI,Trial_I) = max(GROUP2_T{MouseI,Trial_I}.Grid_module);
            Num_Gp2_M1(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 1));
            Num_Gp2_M2(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 2));
            Num_Gp2_M3(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 3));
        catch
        end
    end
end

%% Grid Orientation
%Group2
GridModNaoki =  load('../08_GridModule/GridMod.mat');
A = GridModNaoki ;
ind = ~isnan(A.gp2_GSratio);
Gp2Ori = A.gp2_GSpeak(ind);
Gp2Ori = cell2mat(cellfun(@(x) x(:,2)', Gp2Ori, 'UniformOutput', false));
Group2Diff = abs(Gp2Ori(:,1) - Gp2Ori(:,2));

ind = ~isnan(A.gp1_GSratio);
Gp1Ori = A.gp1_GSpeak(ind);
Gp1Ori = cell2mat(cellfun(@(x) x(:,2)', Gp1Ori, 'UniformOutput', false));
Group1Diff = abs(Gp1Ori(:,1) - Gp1Ori(:,2));

[~, p_Group2] = ttest2(Group2Diff,Group1Diff);


means = [mean(Group1Diff,'omitnan'), mean(Group2Diff,'omitnan')];
sems  = [std(Group1Diff,'omitnan')/sqrt(sum(~isnan(Group1Diff))), ...
    std(Group2Diff,'omitnan')/sqrt(sum(~isnan(Group2Diff)))];

figure;
b = bar(means); hold on;
errorbar(1:2, means, sems, 'k','LineStyle','none','LineWidth',1.5);
set(gca,'XTick',1:2,'XTickLabel',{'Group1','Group2'})
ylabel('Grid Orientation Difference (deg)')

ymax = max(means + sems)*1.1;
line([1 2],[ymax ymax],'Color','k','LineWidth',1.5)
if p_Group2 < 0.001
    text(1.5, ymax, '***','HorizontalAlignment','center','FontSize',14)
elseif p_Group2 < 0.01
    text(1.5, ymax, '**','HorizontalAlignment','center','FontSize',14)
elseif p_Group2 < 0.05
    text(1.5, ymax, '*','HorizontalAlignment','center','FontSize',14)
else
    text(1.5, ymax, 'n.s.','HorizontalAlignment','center')
end

title('Grid Orientation Difference: Group2 vs Group1')
box off

%%
A = GridModNaoki;


%% ========== Compute Results of grid orientation ==========

% ref = true;
ref = false;

Results_Group1  = compute_results_GridOri(GROUP1_T, A.gp1_GSpeak,  min_modulenum, mode, ref);
Results_Group2 = compute_results_GridOri(GROUP2_T, A.gp2_GSpeak, min_modulenum, mode, ref);

% ========== Prepare metrics ==========
metrics = {'within_mean','between_mean','within_std','between_std','within_IQR','between_IQR'};
Data_Group1  = extract_metrics(Results_Group1, metrics);
Data_Group2 = extract_metrics(Results_Group2, metrics);

% ========== CV per module ==========
Data_Group1.CV   = Data_Group1.within_std ./ Data_Group1.within_mean;
Data_Group2.CV  = Data_Group2.within_std ./ Data_Group2.within_mean;

% ========== Statistical tests ==========
[p.CV] = compare_groups(Data_Group1.CV, Data_Group2.CV);


%% ===================
within_Group1  = [Results_Group1.within_mean]';
between_Group1 = [Results_Group1.between_mean]';

within_Group2  = [Results_Group2.within_mean]';
between_Group2 = [Results_Group2.between_mean]';

% CV (within-module)
CV_Group1  = [Results_Group1.within_std] ./ [Results_Group1.within_mean];
CV_Group2 = [Results_Group2.within_std] ./ [Results_Group2.within_mean];

CV_Group1  = CV_Group1(~isnan(CV_Group1));
CV_Group2 = CV_Group2(~isnan(CV_Group2));

%% =================== ANOVA: Within/Between ===================
data = [within_Group1; between_Group1; within_Group2; between_Group2];
group = [ ...
    repmat("Group1", length(within_Group1),1);
    repmat("Group1", length(between_Group1),1);
    repmat("Group2", length(within_Group2),1);
    repmat("Group2", length(between_Group2),1)];
condition = [ ...
    repmat("Within", length(within_Group1),1);
    repmat("Between", length(between_Group1),1);
    repmat("Within", length(within_Group2),1);
    repmat("Between", length(between_Group2),1)];

[p_anova, tbl, stats] = anovan(data, {group, condition}, ...
    'model','interaction', ...
    'varnames',{'Group','Condition'});

means = [
    mean(within_Group1,'omitnan'), mean(between_Group1,'omitnan');
    mean(within_Group2,'omitnan'), mean(between_Group2,'omitnan')];

sems = [
    std(within_Group1,'omitnan')/sqrt(sum(~isnan(within_Group1))), ...
    std(between_Group1,'omitnan')/sqrt(sum(~isnan(between_Group1)));
    std(within_Group2,'omitnan')/sqrt(sum(~isnan(within_Group2))), ...
    std(between_Group2,'omitnan')/sqrt(sum(~isnan(between_Group2)))];

%% =================== Within/Between ===================
figure;
b = bar(means); hold on;

[ngroups, nbars] = size(means);
for i = 1:nbars
    x = (1:ngroups) - 0.15 + (i-1)*0.3;
    errorbar(x, means(:,i), sems(:,i), 'k','LineStyle','none','LineWidth',1.5);
end

set(gca,'XTick',1:2,'XTickLabel',{'Group1','Group2'})
legend({'Within','Between'})
ylabel('Grid orientation difference (deg)')
title('Grid orientation difference')
box off

% ===================
figure
c = multcompare(stats, 'Dimension',[1 2]);

gnames = {'Group1-Within','Group1-Between','Group2-Within','Group2-Between'};
for i = 1:size(c,1)
    g1 = gnames{c(i,1)};
    g2 = gnames{c(i,2)};
    fprintf('%s vs %s : diff=%.2f, p=%.3f\n', g1, g2, c(i,3), c(i,6));
end


%% Fig

CGroup1 = [1 1 1];
CGroup2 = [255 168 60]/255;

pos = [300 500 130 160];
fig = figure('Position', pos);


Group1 = within_Group1;
Group2 = within_Group2;

means = [
    mean(Group1,'omitnan')
    mean(Group2,'omitnan')
    ];

sems = [
    std(Group1,'omitnan')/sqrt(sum(~isnan(Group1)))
    std(Group2,'omitnan')/sqrt(sum(~isnan(Group2)))
    ];

x = [0.9, 2.1];
b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
hold on

b.CData(1,:) = [1 1 1];
b.CData(2,:) = CGroup2;

errorbar(x, means, sems, ...
    'k','LineStyle','none','LineWidth',0.5,'CapSize',12)


ytol = 0.02;
dx = 0.15;
x1 = SimpleBeeSwarm(Group1, x(1), ytol, dx);
x2 = SimpleBeeSwarm(Group2, x(2), ytol, dx);
scatter(x1, Group1, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
scatter(x2, Group2, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)

% xlim([0.5 2.5])
xlim([0 3])

set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{})

box off



A = Group1(~isnan(Group1));
B = Group2(~isnan(Group2));


[~, p_ttest, ~, stats] = ttest2(A,B);
tval = stats.tstat;
df   = stats.df;
fprintf('p = %.4f, t(%d) = %.3f\n', p_ttest, df, tval);


yticks([0 3 6])
ylim([0 6])
ylabel('')
xticks([])
yticklabels([])


s = strcat("FigS_GridOriDiff_Group2.emf");
exportgraphics(fig, s);

% A = Y{1};
%     B = Y{2};
mA = mean(A);
mB = mean(B);
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
fprintf('Group2: %.3f ± %.3f\n', mB, sdB);

%% =================== Functions ===================
function Data = extract_metrics(Results, metrics)
for i = 1:length(metrics)
    mname = metrics{i};
    Data.(mname) = [Results.(mname)]';
end
end

%%
% function Results = compute_results_GridOri(OcT, OriPeak_all, min_modulenum, mode, use_ref_module)
function Results = compute_results_GridOri(OcT, GSpeak, min_modulenum, mode, use_ref_module)

% Results = struct( ...
%     'within_mean', {}, ...
%     'between_mean', {}, ...
%     'within_std', {}, ...
%     'between_std', {}, ...
%     'within_IQR', {}, ...
%     'between_IQR', {});
Results = struct([]);


for m = 1:size(OcT,1)
    for d = 1:size(OcT,2)

        tmp = OcT{m,d};
        if isempty(tmp) || height(tmp)==0
            continue
        end

        % ===== GSpeakから直接取得 =====
        peak_cell = GSpeak{m,d};
        if isempty(peak_cell) || any(isnan(peak_cell(:)))
            continue
        end
        module_peak = peak_cell(:,2)'; % 例: Nx2 の2列目がorientation

        if length(module_peak) < min_modulenum % module数チェック（安全）
            continue
        end


        % ===== Grid情報 =====
        Gori = tmp.Grid_orientaion;
        Gmod = tmp.Grid_module;

        idx = (Gori~=0) & (Gmod>0);
        if sum(idx) < 1
            continue
        end

        ori = Gori(idx);
        module_id = Gmod(idx);

        unique_mod = unique(module_id);
        if length(unique_mod) < min_modulenum
            continue
        end


        % ===== 計算 =====
        out = compute_orientation_difference_peak(ori, module_id, module_peak, use_ref_module);
        Results = [Results; out];

    end
end
end



%%
function out = compute_orientation_difference_peak(ori, module_id, module_peak, use_ref_module)

% unique_mod = sort(unique(module_id));
% nMod = length(unique_mod);
unique_mod = [1 2];
nMod = 2;

within = [];
between = [];

% ===== within =====
for i = 1:nMod
    idx_i = module_id == unique_mod(i);
    oi = ori(idx_i);

    peak_i = module_peak(i);

    diff = angle_diff(oi, peak_i);
    within = [within; diff(:)];
end

% ===== between =====
if use_ref_module
    ref_peak = module_peak(1);

    for i = 2:nMod
        idx_i = module_id == unique_mod(i);
        oi = ori(idx_i);

        diff = angle_diff(oi, ref_peak);
        between = [between; diff(:)];
    end

else
    for i = 1:nMod
        idx_i = module_id == unique_mod(i);
        oi = ori(idx_i);

        other_peaks = module_peak;
        other_peaks(i) = [];

        for j = 1:length(other_peaks)
            diff = angle_diff(oi, other_peaks(j));
            between = [between; diff(:)];
        end
    end
end

% ===== summary =====
out.within = within;
out.between = between;

out.within_mean  = mean(within,'omitnan');
out.between_mean = mean(between,'omitnan');

out.within_std   = std(within,'omitnan');
out.between_std  = std(between,'omitnan');

out.within_IQR   = iqr(within);
out.between_IQR  = iqr(between);
end


%%
function d = angle_diff(a, b)
d = abs(a - b);
d = mod(d,180);
d(d>90) = 180 - d(d>90);
end


%%
function [p] = compare_groups(data1, data2)
% NaN除外
data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

[~,p] = ttest2(data1,data2);
end


%%
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