%% ========== Load data ==========
clear; close all hidden


GROUP1 = load("../Group1_CellTables.mat");
GROUP2 = load("../Group2_CellTables.mat");
GROUP1_T = GROUP1.GROUP1_T;
GROUP2_T = GROUP2.GROUP2_T;

min_modulenum = 2;
% min_modulenum = 1;
mode = 'adjacent'; % 'all' 


%% Use only module 1 and 2
GROUP1_T_modified = GROUP1_T;
GROUP2_T_modified = GROUP2_T;

for m = 1:numel(GROUP1_T_modified)
    if ~isempty(GROUP1_T_modified{m})
        idx3 = GROUP1_T_modified{m}.Grid_module == 3;
        GROUP1_T_modified{m}.Grid_module(idx3) = 0;
    end
end

for m = 1:numel(GROUP2_T_modified)
    if ~isempty(GROUP2_T_modified{m})
        idx3 = GROUP2_T_modified{m}.Grid_module == 3;
        GROUP2_T_modified{m}.Grid_module(idx3) = 0;
    end
end

GROUP1_T = GROUP1_T_modified;
GROUP2_T = GROUP2_T_modified;


%% max grid module number
MaxModule_Grp1 = nan(5,3);
Num_Grp1_M1 = nan(5,3);
Num_Grp1_M2 = nan(5,3);
Num_Grp1_M3 = nan(5,3);
for MouseI = 1:5
    for Trial_I = 1:3
        MaxModule_Grp1(MouseI,Trial_I) = max(GROUP1_T{MouseI,Trial_I}.Grid_module);
        Num_Grp1_M1(MouseI,Trial_I) = length(find(GROUP1_T{MouseI,Trial_I}.Grid_module == 1));
        Num_Grp1_M2(MouseI,Trial_I) = length(find(GROUP1_T{MouseI,Trial_I}.Grid_module == 2));
        Num_Grp1_M3(MouseI,Trial_I) = length(find(GROUP1_T{MouseI,Trial_I}.Grid_module == 3));
    end
end

MaxModule_Grp2 = nan(7,3);
Num_Grp2_M1 = nan(7,3);
Num_Grp2_M2 = nan(7,3);
Num_Grp2_M3 = nan(7,3);
for MouseI = 1:7
    for Trial_I = 1:3
        try
            MaxModule_Grp2(MouseI,Trial_I) = max(GROUP2_T{MouseI,Trial_I}.Grid_module);
            Num_Grp2_M1(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 1));
            Num_Grp2_M2(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 2));
            Num_Grp2_M3(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 3));
        catch
        end
    end
end

%% GSratio

GridModNaoki =  load("08_GridModule\GridMod.mat");


A = GridModNaoki ;


GP2ratio = A.gp2_GSratio(:);
GP2ratio = GP2ratio(~isnan(GP2ratio));
GP1ratio = A.gp1_GSratio(:);
GP1ratio = GP1ratio(~isnan(GP1ratio));

[~, p_GP2] = ttest2(GP2ratio, GP1ratio);



%% ========== Compute Results of grid scale ==========

scale_mode = 'ratio'; % 
% scale_mode = 'raw'; % not scale ratio but raw scale [cm]
Results_Grp1  = compute_results_GridScale(GROUP1_T, A.gp1_GSpeak, A.gp1_GSratio, min_modulenum, mode, scale_mode);
Results_Grp2 = compute_results_GridScale(GROUP2_T, A.gp2_GSpeak, A.gp2_GSpeak, min_modulenum, mode, scale_mode);



%% ==================
within_Grp1  = [Results_Grp1.within_mean]';
between_Grp1 = [Results_Grp1.between_mean]';

within_Grp2  = [Results_Grp2.within_mean]';
between_Grp2 = [Results_Grp2.between_mean]';



% %% statistics
% =================== ANOVA: Within/Between ===================
data = [within_Grp1; between_Grp1; within_Grp2; between_Grp2];
group = [ ...
    repmat("Group1", length(within_Grp1),1);
    repmat("Group1", length(between_Grp1),1);
    repmat("Group2", length(within_Grp2),1);
    repmat("Group2", length(between_Grp2),1)];
condition = [ ...
    repmat("Within", length(within_Grp1),1);
    repmat("Between", length(between_Grp1),1);
    repmat("Within", length(within_Grp2),1);
    repmat("Between", length(between_Grp2),1)];

[p_anova, tbl, stats] = anovan(data, {group, condition}, ...
    'model','interaction', ...
    'varnames',{'Group','Condition'});

means = [
    mean(within_Grp1,'omitnan'), mean(between_Grp1,'omitnan');
    mean(within_Grp2,'omitnan'), mean(between_Grp2,'omitnan')];

sems = [
    std(within_Grp1,'omitnan')/sqrt(sum(~isnan(within_Grp1))), ...
    std(between_Grp1,'omitnan')/sqrt(sum(~isnan(between_Grp1)));
    std(within_Grp2,'omitnan')/sqrt(sum(~isnan(within_Grp2))), ...
    std(between_Grp2,'omitnan')/sqrt(sum(~isnan(between_Grp2)))];


df_error = tbl{5,3};
F_group = tbl{2,6};
p_group = tbl{2,7};
df_group = tbl{2,3};
F_cond = tbl{3,6};
p_cond = tbl{3,7};
df_cond = tbl{3,3};
F_inter = tbl{4,6};
p_inter = tbl{4,7};
df_inter = tbl{4,3};
fprintf('Group × Condition: F(%d,%d)=%.3f, p=%.4f\n', ...
    df_inter, df_error, F_inter, p_inter)
fprintf('Group: F(%d,%d)=%.3f, p=%.4f\n', ...
    df_group, df_error, F_group, p_group)
fprintf('Condition: F(%d,%d)=%.3f, p=%.4g\n', ...
    df_cond, df_error, F_cond, p_cond)



%% plot adjusting bar-spacing (improved)

% Cgroup1 = [1 1 1];
% CGroup2 = [255 168 60]/255;
Cgroup1 = [1 1 1];
CGroup2 = [240 134 134]/255;

% close all

pos = [100 700 150 111];
fig = figure;
set(fig, 'Position', pos)

hold on


groupGap  = 0.45;   % ← Group1 vs Group2 
withinGap = 0.18;   % ← Within vs Between 
barWidth  = 0.14;   % 

% =========================================================
% x positions
% =========================================================

ctrlCenter = 1;
Group2Center = ctrlCenter + groupGap;

x_ctrl_within  = ctrlCenter - withinGap/2;
x_ctrl_between = ctrlCenter + withinGap/2;

x_Grp2_within  = Group2Center - withinGap/2;
x_Grp2_between = Group2Center + withinGap/2;

% =========================================================
% bar
% =========================================================

bar(x_ctrl_within, means(1,1), ...
    barWidth, ...
    'FaceColor','w', ...
    'EdgeColor','k')

bar(x_ctrl_between, means(1,2), ...
    barWidth, ...
    'FaceColor','w', ...
    'EdgeColor','k')

bar(x_Grp2_within, means(2,1), ...
    barWidth, ...
    'FaceColor',CGroup2, ...
    'EdgeColor','k')

bar(x_Grp2_between, means(2,2), ...
    barWidth, ...
    'FaceColor',CGroup2, ...
    'EdgeColor','k')

% =========================================================
% error bar
% =========================================================

errorbar(x_ctrl_within, means(1,1), sems(1,1), ...
    'k', 'LineStyle','none', ...
    'LineWidth',0.5, 'CapSize',8)

errorbar(x_ctrl_between, means(1,2), sems(1,2), ...
    'k', 'LineStyle','none', ...
    'LineWidth',0.5, 'CapSize',8)

errorbar(x_Grp2_within, means(2,1), sems(2,1), ...
    'k', 'LineStyle','none', ...
    'LineWidth',0.5, 'CapSize',8)

errorbar(x_Grp2_between, means(2,2), sems(2,2), ...
    'k', 'LineStyle','none', ...
    'LineWidth',0.5, 'CapSize',8)

% =========================================================
% scatter
% =========================================================

% ytol = 0.03;
ytol = 0.01;
dx = 0.02;

% --- Group1 ---
x1 = SimpleBeeSwarm(within_Grp1, x_ctrl_within, ytol, dx);
x2 = SimpleBeeSwarm(between_Grp1, x_ctrl_between, ytol, dx);

scatter(x1, within_Grp1, 8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',0.7)

scatter(x2, between_Grp1, 8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',0.7)

% --- Group2 ---
x1 = SimpleBeeSwarm(within_Grp2, x_Grp2_within, ytol, dx);
x2 = SimpleBeeSwarm(between_Grp2, x_Grp2_between, ytol, dx);

scatter(x1, within_Grp2, 8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',0.7)

scatter(x2, between_Grp2, 8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','w', ...
    'LineWidth',0.7)

% =========================================================
% axis
% =========================================================

xlim([0.75 Group2Center+0.25])

ylim([0 1])

yticks([0 .5 1])
yticklabels([])

set(gca,'XTick',[])
box off

xlabel('')
ylabel('')

set(gca,'LineWidth',1)


A = within_Grp1;
B = between_Grp1;
mA = mean(A);
mB = mean(B);
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
fprintf('Group2: %.3f ± %.3f\n', mB, sdB);

A = within_Grp2;
B = between_Grp2;
mA = mean(A);
mB = mean(B);
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
fprintf('Group2: %.3f ± %.3f\n', mB, sdB);


s = strcat("GridScale_2.emf");
exportgraphics(fig, s);



[~, p_GP2_within] = ttest2(within_Grp1, within_Grp2);

%%
figure
c = multcompare(stats, 'Dimension',[1 2]);

% gnames = {'Group1-Within','Group1-Between','Group2-Within','Group2-Between'};r
gnames = {'Group1-Within','Group1-Between','Group2-Within','Group2-Between'};
for i = 1:size(c,1)
    g1 = gnames{c(i,1)};
    g2 = gnames{c(i,2)};
    fprintf('%s vs %s : diff=%.2f, p=%.3f\n', g1, g2, c(i,3), c(i,6));
end

%% % q value
figure
% c = multcompare(stats, 'Dimension',[1 2]);
% gnames = {'Group1-Within','Group1-Between', ...
%           'Group2-Within','Group2-Between'};

% Error term from ANOVA table
MSE = tbl{5,5};      % Mean Sq. of Error
df_error = tbl{5,3}; % Error d.f.

% sample sizes
n = [
    length(within_Grp1)
    length(between_Grp1)
    length(within_Grp2)
    length(between_Grp2)
];

for i = 1:size(c,1)

    g1i = c(i,1);
    g2i = c(i,2);

    g1 = gnames{g1i};
    g2 = gnames{g2i};

    diffmean = abs(c(i,4));
    p = c(i,6);

    % harmonic mean for unequal n
    nh = 2/(1/n(g1i) + 1/n(g2i));

    % Tukey q statistic
    q = diffmean / sqrt(MSE/nh);

    % p formatting
    if p < 0.0001
        ptxt = 'p < 0.0001';
    else
        ptxt = sprintf('p = %.4f', p);
    end

    fprintf('%s vs. %s: %s, q(%d) = %.3f\n', ...
        g1, g2, ptxt, df_error, q);

end


%%
%% =========================================================
% One-way ANOVA
% 4 groups:
% Group1-Within
% Group1-Between
% Group2-Within
% Group2-Between
%% =========================================================

% data
data = [ ...
    within_Grp1;
    between_Grp1;
    within_Grp2;
    between_Grp2];

% 4-level group label
group4 = [ ...
    repmat("Group1-Within",    length(within_Grp1),1);
    repmat("Group1-Between",   length(between_Grp1),1);
    repmat("Group2-Within",   length(within_Grp2),1);
    repmat("Group2-Between",  length(between_Grp2),1)];

%% =========================================================
% One-way ANOVA
% =========================================================

[p_anova, tbl, stats] = anova1(data, group4, 'off');

% =========================================================
% Tukey post hoc
% =========================================================

results = multcompare(stats, ...
    'CType','tukey-kramer', ...
    'Display','off');

% =========================================================
% Mean and SEM
% =========================================================

means = [ ...
    mean(within_Grp1,'omitnan');
    mean(between_Grp1,'omitnan');
    mean(within_Grp2,'omitnan');
    mean(between_Grp2,'omitnan')];

sems = [ ...
    std(within_Grp1,'omitnan')/sqrt(sum(~isnan(within_Grp1)));
    std(between_Grp1,'omitnan')/sqrt(sum(~isnan(between_Grp1)));
    std(within_Grp2,'omitnan')/sqrt(sum(~isnan(within_Grp2)));
    std(between_Grp2,'omitnan')/sqrt(sum(~isnan(between_Grp2)))];

% =========================================================
% ANOVA statistics
% =========================================================

df_group = tbl{2,3};
df_error = tbl{3,3};

F_value = tbl{2,5};
p_value = tbl{2,6};

fprintf('\n====================================\n')
fprintf('One-way ANOVA\n')
fprintf('====================================\n')

fprintf('F(%d,%d)=%.3f, p=%.4g\n', ...
    df_group, df_error, F_value, p_value)

% =========================================================
% Tukey results
% =========================================================

group_names = stats.gnames;

fprintf('\n====================================\n')
fprintf('Tukey multiple comparisons\n')
fprintf('====================================\n')

for i = 1:size(results,1)

    g1 = results(i,1);
    g2 = results(i,2);

    p = results(i,6);

    % Tukey q statistic
    q = abs(results(i,4)) / results(i,5);

    fprintf('%s vs. %s: p = %.4f, q(%d) = %.3f\n', ...
        group_names{g1}, ...
        group_names{g2}, ...
        p, ...
        df_error, ...
        q)

end




%% functions
%% compute grid scale ratio
function Results = compute_results_GridScale(OcT, GSpeak_all, GSratio_all, min_modulenum, mode, scale_mode)
%%
% Results = [];
Results = struct([]);

for m = 1:size(OcT,1)
    for d = 1:size(OcT,2)
        % disp(m)
        % disp(d)
        tmp = OcT{m,d};
        if isempty(tmp) || height(tmp)==0
            continue
        end

        % [scale_ratio, module_peak_scale] = func_grid_scale(OcT,m,d);
        % [scale_ratio, module_peak_scale] = func_grid_scale(OcT,m,d);
        Gscale = OcT{m,d}.Grid_scale;
        Gmod   = OcT{m,d}.Grid_module;
        idx = (Gscale>0) & (Gmod>0);

        % ===== GSpeak取得 =====
        GSpeak = GSpeak_all{m,d};
        if isempty(GSpeak) || isnan(GSratio_all(m,d))
            continue
        end
        module_peak_scale = GSpeak(:,1);
        if length(module_peak_scale) < min_modulenum
            continue
        end


        % ===== scale 切り替え =====
        switch scale_mode
            case 'raw'
                scale = Gscale(idx);
            case 'ratio'
                [scale_ratio, DV, MI] = func_grid_scale(OcT, m, d, module_peak_scale); % ← ratioだけ返す
                scale = scale_ratio;
            otherwise
                error('Unknown scale_mode')
        end
        mod = Gmod(idx);
        out = compute_scale_difference(scale, mod, module_peak_scale, mode); % between-cell differences


        if length(module_peak_scale) >= 2 && all(~isnan(module_peak_scale(1:2)))
            out.peaks = (module_peak_scale(1:2) ./ min(module_peak_scale(1:2)))';
        else
            out.peaks = [NaN NaN];
        end


        if isempty(out)
            continue
        end

        out.mouse = m;
        out.day = d;

        out.scale_ratio = scale;
        out.DV = DV;
        out.MI = MI;
        % out.DVnorm = DV - mean(DV); % chenge normalize by mean of each module


        Results = [Results; out];
        % if ~isempty(out)
        %     Results(end+1) = out;
        % end
    end
end

A = vertcat(Results.peaks);
mean(vertcat(Results.peaks));
end



%% Grid scale ratio
function [scale_ratio_global, DV, MI] = func_grid_scale(OcT, Mind, Day, module_peak_scale)

Gscale = OcT{Mind,Day}.Grid_scale;
Gmod   = OcT{Mind,Day}.Grid_module;

idx = (Gscale > 0) & (Gmod > 0);

scale = Gscale(idx);
% mod   = Gmod(idx);

% % moduleごとのpeak（medianでOK：ratio用なので）
% unique_mod = unique(mod);
% module_peak_scale = nan(length(unique_mod),1);
%
% for i = 1:length(unique_mod)
%     module_peak_scale(i) = median(scale(mod == unique_mod(i)));
% end

if length(module_peak_scale) < 2
    scale_ratio_global = [];
    return
end

min_scale = min(module_peak_scale);
scale_ratio_global = scale / min_scale;


DV = OcT{Mind,Day}.DV_position_um(idx);
MI = OcT{Mind,Day}.Grid_module(idx);
end

%%
function out = compute_scale_difference(scale, mod, module_peak_scale, mode)

% mode:
% 'all' 
% 'adjacent' 

unique_mod = unique(mod);
nMod = length(unique_mod);

if strcmp(mode,'adjacent') && nMod >= 2
    [~, order] = sort(module_peak_scale);
    selected_mod = unique_mod(order(1:2));
else
    selected_mod = unique_mod;
end

within = [];
between = [];

for i = 1:length(selected_mod)
    idx_i = mod == selected_mod(i);
    si = scale(idx_i);

    % --- within ---
    if length(si) > 1
        D = abs(si - si');
        within = [within; D(triu(true(size(D)),1))];
    end

    % --- between ---
    for j = i+1:length(selected_mod)
        idx_j = mod == selected_mod(j);
        sj = scale(idx_j);

        [I,J] = meshgrid(si, sj);
        D = abs(I - J);
        between = [between; D(:)];
    end
end

if isempty(within) || isempty(between)
    out = [];
    return
end

% --- summary ---
out.within_mean = mean(within,'omitnan');
out.between_mean = mean(between,'omitnan');

out.within_std = std(within,'omitnan');
out.between_std = std(between,'omitnan');

out.within_IQR = iqr(within);
out.between_IQR = iqr(between);

% out.peaks = [];

end

%%
function [field_raw, field_norm, module_peak_field] = func_grid_field(OcT, Mind, Day)

Gfield = OcT{Mind,Day}.Grid_field;
Gmod   = OcT{Mind,Day}.Grid_module;

idx = (Gfield > 0) & (Gmod > 0);

field = Gfield(idx);
mod   = Gmod(idx);

% unique_mod = unique(mod);
unique_mod = [1 2];
nMod = length(unique_mod);

module_peak_field = nan(nMod,1);

% ===== module peak =====
for i = 1:nMod
    m = unique_mod(i);
    idx_m = (mod == m);
    % module_peak_field(i) = median(field(idx_m));
    module_peak_field(i) = mean(field(idx_m));
end

% ===== raw =====
field_raw = field;


field_norm = field ./ module_peak_field(1);

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