%% ========== Load data ==========
clear; close all hidden

GROUP1 = load("../Group1_CellTables.mat");
GROUP2 = load("../Group2_CellTables.mat");
GROUP1_T = GROUP1.GROUP1_T;
GROUP2_T = GROUP2.GROUP2_T;


Gphase = load("../12_G_PhaseDiff/Grid_phase.mat");

Gphase_Group1  = Gphase.AllP_GridPhase_GROUP1;
Gphase_Group2 = Gphase.AllP_GridPhase_GROUP2;

s = {'CellID_1', 'CellID_2', 'PhaseDistance', 'PhysicalDistance', 'ModuleInd_Cell1', 'ModuleInd_Cell2', 'WithinOrBetween'};


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

MaxModule_Group1 = nan(7,3);
Num_Group1_M1 = nan(7,3);
Num_Group1_M2 = nan(7,3);
Num_Group1_M3 = nan(7,3);
for MouseI = 1:7
    for Trial_I = 1:3
        try
            MaxModule_Group1(MouseI,Trial_I) = max(GROUP2_T{MouseI,Trial_I}.Grid_module);
            Num_Group1_M1(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 1));
            Num_Group1_M2(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 2));
            Num_Group1_M3(MouseI,Trial_I) = length(find(GROUP2_T{MouseI,Trial_I}.Grid_module == 3));
        catch
        end
    end
end

%% =================== Table ===================
[mouseN, DayN] = size(GROUP1_T);
Gphase_Group1_table = cell(mouseN, DayN);
for m = 1:mouseN
    for d = 1:DayN
        try
            CellID_1 = Gphase_Group1{m, d}(:,1);
            CellID_2 = Gphase_Group1{m, d}(:,2);

            [tf1, loc1] = ismember(CellID_1, GROUP1_T{m,d}.Cell_ID);
            [tf2, loc2] = ismember(CellID_2, GROUP1_T{m,d}.Cell_ID);
            M1 = nan(size(CellID_1));
            M1(tf1) = GROUP1_T{m,d}.Grid_module(loc1(tf1));
            M2 = nan(size(CellID_2));
            M2(tf2) = GROUP1_T{m,d}.Grid_module(loc2(tf2));

            % Lg = M1 == M2;
            Lg = (M1 == M2) & ~isnan(M1) & ~isnan(M2);

            Gphase_Group1_table{m, d} = table(CellID_1, CellID_2, Gphase_Group1{m, d}(:,3), Gphase_Group1{m, d}(:,5), M1, M2, Lg,...
                'VariableNames', s);
        catch ME
            warning('Error at m=%d, d=%d: %s', m, d, ME.message);
        end
    end
end

[mouseN, DayN] = size(GROUP2_T);
Gphase_Group2_table = cell(mouseN, DayN);
for m = 1:mouseN
    for d = 1:DayN
        try
            CellID_1 = Gphase_Group2{m, d}(:,1);
            CellID_2 = Gphase_Group2{m, d}(:,2);


            [tf1, loc1] = ismember(CellID_1, GROUP2_T{m,d}.Cell_ID);
            [tf2, loc2] = ismember(CellID_2, GROUP2_T{m,d}.Cell_ID);
            M1 = nan(size(CellID_1));
            M1(tf1) = GROUP2_T{m,d}.Grid_module(loc1(tf1));
            M2 = nan(size(CellID_2));
            M2(tf2) = GROUP2_T{m,d}.Grid_module(loc2(tf2));

            % Lg = M1 == M2;
            Lg = (M1 == M2) & ~isnan(M1) & ~isnan(M2);

            Gphase_Group2_table{m, d} = table(CellID_1, CellID_2, Gphase_Group2{m, d}(:,3), Gphase_Group2{m, d}(:,5), M1, M2, Lg,...
                'VariableNames', s);
        catch ME
            warning('Error at m=%d, d=%d: %s', m, d, ME.message);
        end
    end
end

%% =================== Compute (within only) ===================
Results_Group1  = compute_within_between_phase(Gphase_Group1_table);
Results_Group2 = compute_within_between_phase(Gphase_Group2_table);

%% =================== Phase variability ===================
CV_Group1  = [Results_Group1.within_std] ./ [Results_Group1.within_mean];
CV_Group2 = [Results_Group2.within_std] ./ [Results_Group2.within_mean];

CV_Group1  = CV_Group1(~isnan(CV_Group1));
CV_Group2 = CV_Group2(~isnan(CV_Group2));

[p_cv] = compare_groups(CV_Group1, CV_Group2);

figure
plot_bar_with_sem({CV_Group1, CV_Group2}, p_cv, {'Group1','Group2'}, 'CV of within-module phase');


%% Between module analysis  plot
within_Group1  = [Results_Group1.within_mean];
between_Group1 = [Results_Group1.between_mean];
within_Group2  = [Results_Group2.within_mean];
between_Group2 = [Results_Group2.between_mean];

data = [within_Group1, between_Group1, within_Group2, between_Group2]';
Group = [ ...
    repmat({'Group1'}, length(within_Group1),1);
    repmat({'Group1'}, length(between_Group1),1);
    repmat({'Group2'},    length(within_Group2),1);
    repmat({'Group2'},    length(between_Group2),1)];
WB = [ ...
    repmat({'Within'},  length(within_Group1),1);
    repmat({'Between'}, length(between_Group1),1);
    repmat({'Within'},  length(within_Group2),1);
    repmat({'Between'}, length(between_Group2),1)];

[p,tbl,stats] = anovan(data, {Group, WB}, ...
    'model','interaction', ...
    'varnames',{'Group','WithinBetween'}); % see Group × WithinBetween
figure;
[c,~,~,gnames] = multcompare(stats, 'Dimension',[1 2], 'CType','tukey-kramer');

%% plot
means_raw = [
    mean(within_Group1,'omitnan'), mean(between_Group1,'omitnan');
    mean(within_Group2,'omitnan'), mean(between_Group2,'omitnan')];

sems_raw = [
    std(within_Group1,'omitnan')/sqrt(sum(~isnan(within_Group1))), ...
    std(between_Group1,'omitnan')/sqrt(sum(~isnan(between_Group1)));
    std(within_Group2,'omitnan')/sqrt(sum(~isnan(within_Group2))), ...
    std(between_Group2,'omitnan')/sqrt(sum(~isnan(between_Group2)))];

figure; hold on;

b = bar(means_raw);

[ngroups, nbars] = size(means_raw);
for i = 1:nbars
    x = (1:ngroups) - 0.15 + (i-1)*0.3;
    errorbar(x, means_raw(:,i), sems_raw(:,i), ...
        'k','LineStyle','none','LineWidth',1.5);
end

set(gca,'XTick',1:2,'XTickLabel',{'Group1','Group2'})
% legend({'Within','Between'})
ylabel('Phase distance')
title('Within vs Between (Grid phase)')
box off




%% =================== Distance–Phase (Binning) ===================
bin_width = 10; % um
max_dist = 500;
edges = 0:bin_width:max_dist;
bin_centers = edges(1:end-1) + bin_width/2;

% --- Compute binned phase for each type ---
[mean_Group1_within, sem_Group1_within, allbins_Group1_within]   = compute_binned_phase_module(Results_Group1, edges, 'within');
[mean_Group1_between, sem_Group1_between, allbins_Group1_between] = compute_binned_phase_module(Results_Group1, edges, 'between');

[mean_Group2_within, sem_Group2_within, allbins_Group2_within]   = compute_binned_phase_module(Results_Group2, edges, 'within');
[mean_Group2_between, sem_Group2_between, allbins_Group2_between] = compute_binned_phase_module(Results_Group2, edges, 'between');

mean_Group1_within   = interp1(bin_centers(~isnan(mean_Group1_within)), mean_Group1_within(~isnan(mean_Group1_within)), bin_centers, 'linear', NaN);
mean_Group1_between  = interp1(bin_centers(~isnan(mean_Group1_between)), mean_Group1_between(~isnan(mean_Group1_between)), bin_centers, 'linear', NaN);
mean_Group2_within  = interp1(bin_centers(~isnan(mean_Group2_within)), mean_Group2_within(~isnan(mean_Group2_within)), bin_centers, 'linear', NaN);
mean_Group2_between = interp1(bin_centers(~isnan(mean_Group2_between)), mean_Group2_between(~isnan(mean_Group2_between)), bin_centers, 'linear', NaN);


%% =================== Plot ===================
%% within
figure; hold on;

% Group1
h0 = errorbar(bin_centers, mean_Group1_within,  sem_Group1_within,  '-o', 'LineWidth',1, 'Color',[0 0.4 1]);

% Group2
h1 = errorbar(bin_centers, mean_Group2_within,  sem_Group2_within,  '-o', 'LineWidth',1, 'Color',[1 0.2 0.2]);

xlabel('Physical distance (um)')
ylabel('Phase distance')

title('Distance–Phase relationship')
box off
% legend({'Within','Between'})
xlim([0 250])
type = 'within'; %

nShuffle = 100;

shuf_avg = cell(2,1);
shuf_sem = cell(2,1);
shuf_all = cell(2,1);
h_s = cell(2,1);
rng(41)
for G = 1:2
    if G == 1
        R = Results_Group1;
        Cl = [0.6 0.8 1];
    else
        R = Results_Group2;
        Cl = [1 0.6 0.6];
    end

    mean_shuf_mat = nan(length(bin_centers), nShuffle, length(R));  % [bin x shuffle x session]
    nBin = length(edges)-1;
    for i = 1:length(R)
        if strcmp(type,'within')
            phys  = R(i).within_phys;
            phase = R(i).within_phase;
        else
            phys  = R(i).between_phys;
            phase = R(i).between_phase;
        end

        for s = 1:nShuffle
            shuffle_phase = phase(randperm(length(phase)));
            [mean_shuf, ~, all_bins] = compute_binned_phase_single(phys, shuffle_phase, edges);
            mean_shuf_mat(:,s,i) = mean_shuf;  % セッションごと
        end
    end

    result = reshape(mean_shuf_mat, nBin, []);



    shuf_all{G} = result;
    shuf_avg{G} = mean(shuf_all{G},2,'omitnan'); %
    shuf_sem{G} = std(shuf_all{G},0,2,'omitnan')./sqrt(sum(~isnan(shuf_all{G}), 2));

    h_s{G} = errorbar(bin_centers, shuf_avg{G}, shuf_sem{G}, '-o', 'Color', Cl, 'LineWidth',1);

end
% legend({'Group1-shuffle','Group1-Between','Group1-Within','Group2-Within'}, 'Location','best')
legend([h0, h1, h_s{1}, h_s{2}], 'Group1','Group2','Group1-shuffle','Group2-shuffle');
ylim([0 0.6])


%%
p_bin_Group1 = nan(length(bin_centers),1);
for b = 1:length(bin_centers)
    G = 1;
    data_real = allbins_Group1_within{b};
    data_shuf = shuf_all{G}(b,:);
    if ~isempty(data_real) && ~isempty(data_shuf)
        [~,p] = ttest2(data_real, data_shuf); % Welch t-test
        p_bin_Group1(b) = p;
    end
end
p_fdr_Group1 = mafdr(p_bin_Group1,'BHFDR',true);

p_bin_Group1 = nan(length(bin_centers),1);
for b = 1:length(bin_centers)
    G = 2;
    data_real = allbins_Group1_within{b};
    data_shuf = shuf_all{G}(b,:);
    if ~isempty(data_real) && ~isempty(data_shuf)
        [~,p] = ttest2(data_real, data_shuf); % Welch t-test
        p_bin_Group1(b) = p;
    end
end
p_fdr_Group1 = mafdr(p_bin_Group1,'BHFDR',true);

%%

[p_Group1_vs_Group2_fdr, p_Group1_vs_shuf_fdr, p_Group2_vs_shuf_fdr] = ...
    func_fourttest(bin_centers, allbins_Group1_within, allbins_Group2_within, shuf_all);



%% bar plot

i = 2;
c1 = func_phasedistance_bar(i, allbins_Group1_within, allbins_Group2_within, shuf_all);

i = 16;
c2 = func_phasedistance_bar(i, allbins_Group1_within, allbins_Group2_within, shuf_all);




%% =================== Plot ===================
%% between
figure; hold on;

% Group1
h0 = errorbar(bin_centers, mean_Group1_between, sem_Group1_between, '-o', 'LineWidth',1, 'Color',[0 0.4 1]);

% Group2
h1 = errorbar(bin_centers, mean_Group2_between, sem_Group2_between, '-o', 'LineWidth',1, 'Color',[1 0.2 0.2]);

xlabel('Physical distance (um)')
ylabel('Phase distance')
title('Distance–Phase relationship')
box off
xlim([0 250])
%
% type = 'within'; 
type = 'between';


% nShuffle = 100;
mean_shuf_avg = cell(2,1);
mean_shuf_sem = cell(2,1);
mean_shuf_all = cell(2,1);
h_s = cell(2,1);
for G = 1:2
    if G == 1
        R = Results_Group1;
        Cl = [0.6 0.8 1];
    else
        R = Results_Group2;
        Cl = [1 0.6 0.6];
    end

    mean_shuf_mat = nan(length(bin_centers), nShuffle, length(R));  % [bin x shuffle x session]

    for i = 1:length(R)
        if strcmp(type,'within')
            phys  = R(i).within_phys;
            phase = R(i).within_phase;
        else
            phys  = R(i).between_phys;
            phase = R(i).between_phase;
        end

        for s = 1:nShuffle
            shuffle_phase = phase(randperm(length(phase)));
            [mean_shuf, ~] = compute_binned_phase_single(phys, shuffle_phase, edges);
            mean_shuf_mat(:,s,i) = mean_shuf; 
        end
    end


    mean_shuf_all{G} = mean(mean_shuf_mat,2,'omitnan');
    mean_shuf_all{G} = squeeze(mean_shuf_all{G});

    mean_shuf_avg{G} = mean(mean_shuf_all{G},2,'omitnan');
    mean_shuf_sem{G} = std(mean_shuf_all{G},0,2,'omitnan')./sqrt(sum(~isnan(mean_shuf_all{G}), 2));

    h_s{G} = errorbar(bin_centers, mean_shuf_avg{G}, mean_shuf_sem{G}, '-o', 'Color', Cl, 'LineWidth',1);
    uistack(h_s{G},'bottom');

end
legend([h0, h1, h_s{1}, h_s{2}], 'Group1','Group2','Group1-shuffle','Group2-shuffle');
ylim([0 0.6])

%%
p_bin_Group1_between = nan(length(bin_centers),1);
for b = 1:length(bin_centers)
    G = 1;
    data_real = allbins_Group1_within{b};
    data_shuf = mean_shuf_all{G}(b,:);
    if ~isempty(data_real) && ~isempty(data_shuf)
        [~,p] = ttest2(data_real, data_shuf); % Welch t-test
        p_bin_Group1_between(b) = p;
    end
end
p_fdr_Group1_between = mafdr(p_bin_Group1_between,'BHFDR',true);

p_bin_Group1_between = nan(length(bin_centers),1);
for b = 1:length(bin_centers)
    G = 2;
    data_real = allbins_Group1_within{b};
    data_shuf = mean_shuf_all{G}(b,:);
    if ~isempty(data_real) && ~isempty(data_shuf)
        [~,p] = ttest2(data_real, data_shuf); % Welch t-test
        p_bin_Group1_between(b) = p;
    end
end
p_fdr_Group1_between = mafdr(p_bin_Group1_between,'BHFDR',true);



%%
%% =================== All-session scatter (within + between) ===================

all_phys_Group1_within  = [];
all_phase_Group1_within = [];
all_phys_Group1_between = [];
all_phase_Group1_between = [];

for i = 1:length(Results_Group1)
    all_phys_Group1_within  = [all_phys_Group1_within; Results_Group1(i).within_phys(:)];
    all_phase_Group1_within = [all_phase_Group1_within; Results_Group1(i).within_phase(:)];

    all_phys_Group1_between  = [all_phys_Group1_between; Results_Group1(i).between_phys(:)];
    all_phase_Group1_between = [all_phase_Group1_between; Results_Group1(i).between_phase(:)];
end

all_phys_Group2_within  = [];
all_phase_Group2_within = [];
all_phys_Group2_between = [];
all_phase_Group2_between = [];

for i = 1:length(Results_Group2)
    all_phys_Group2_within  = [all_phys_Group2_within; Results_Group2(i).within_phys(:)];
    all_phase_Group2_within = [all_phase_Group2_within; Results_Group2(i).within_phase(:)];

    all_phys_Group2_between  = [all_phys_Group2_between; Results_Group2(i).between_phys(:)];
    all_phase_Group2_between = [all_phase_Group2_between; Results_Group2(i).between_phase(:)];
end

% =================== Scatter plot ===================
figure; hold on;

scatter(all_phys_Group1_within,  all_phase_Group1_within,  5, [0 0.4 1], 'filled', 'MarkerFaceAlpha',0.6)
scatter(all_phys_Group1_between, all_phase_Group1_between, 5, [0.6 0.8 1], 'filled', 'MarkerFaceAlpha',0.3)

scatter(all_phys_Group2_within,  all_phase_Group2_within,  5, [1 0.2 0.2], 'filled', 'MarkerFaceAlpha',0.6)
scatter(all_phys_Group2_between, all_phase_Group2_between, 5, [1 0.6 0.6], 'filled', 'MarkerFaceAlpha',0.3)

xlabel('Physical distance (um)')
ylabel('Phase distance')
legend({'Group1-Within','Group1-Between','Group2-Within','Group2-Between'}, 'Location','best')
title('All sessions pooled (within vs between)')
box off



%% =================== Within-session scatter ===================
Ccont = [0 0 0];
CGroup2 = [255 168 60]/255;

all_phys_Group1_within  = [];
all_phase_Group1_within = [];
all_phys_Group1_between = [];
all_phase_Group1_between = [];

for i = 1:length(Results_Group1)
    all_phys_Group1_within  = [all_phys_Group1_within; Results_Group1(i).within_phys(:)];
    all_phase_Group1_within = [all_phase_Group1_within; Results_Group1(i).within_phase(:)];

end

all_phys_Group2_within  = [];
all_phase_Group2_within = [];


for i = 1:length(Results_Group2)
    all_phys_Group2_within  = [all_phys_Group2_within; Results_Group2(i).within_phys(:)];
    all_phase_Group2_within = [all_phase_Group2_within; Results_Group2(i).within_phase(:)];

end

% =================== Scatter plot ===================
figure; hold on;

% ===== Group1 =====
scatter(all_phys_Group1_within,  all_phase_Group1_within, 8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.5)

% ===== Group2 =====
scatter(all_phys_Group2_within,  all_phase_Group2_within, 8, ...
    'MarkerEdgeColor',CGroup2, ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.5)

xlabel('Physical distance (um)')
ylabel('Phase distance')
% legend({'Group1-Within','Group1-Between','Group2-Within','Group2-Between'}, 'Location','best')
title('All sessions pooled (within)')
box off


%%
%% =================== Violin plot: Phase distance only ===================

all_phase_Group1 = all_phase_Group1_within;
all_phase_Group2 = all_phase_Group2_within;

ydata = [all_phase_Group1; all_phase_Group2];

group = [ ...
    repmat("Group1", length(all_phase_Group1), 1); ...
    repmat("Group2", length(all_phase_Group2), 1)];

% categorical
group = categorical(group);


valid = ~isnan(ydata);
ydata = ydata(valid);
group = group(valid);

% ===== Violin plot =====
figure;
violinplot(group, ydata);

ylabel('Grid phase distance')
title('Phase distance distribution (Group1 vs Group2)')
box off


%% =================== Violin plot (Group1 left, Group2 right) ===================

pos = [300 500 130 160];
fig2 = figure('Position', pos);

Ccont = [0 0 0];
CGroup2 = [255 168 60]/255;

% =====
all_phase_Group1 = all_phase_Group1_within;
all_phase_Group2 = all_phase_Group2_within;

ydata = [all_phase_Group1; all_phase_Group2];





%%

pos = [300 500 130 160];
fig2 = figure('Position', pos);


Y = {all_phase_Group1, all_phase_Group2};

% plot([0 4], [1 1], ':', 'Color', [1 1 1]*0.5)
hold on

xgroupdata = categorical(repelem(["group1";"group2"],[length(Y{1}), length(Y{2})]));
v = violinplot_2(cell2mat(Y'), xgroupdata, 'ShowData', false, 'ViolinColor', [1 1 1]*0.8,...
    'BoxColor', [1 1 1]*0.4, 'Width', 0.3, 'EdgeColor', [1 1 1]*1, 'ViolinAlpha', 1);

v(2).ViolinColor = {CGroup2};

ax = gca;
xlim([0.4 2.6])
box off

xticks('')

ylim([0 1])
yticks([0 0.5 1])
yticklabels({})


[p,h,stats] = ranksum(Y{1}, Y{2});
n1 = numel(Y{1});
n2 = numel(Y{2});
% Mann–Whitney U
U1 = stats.ranksum - n1*(n1+1)/2;
U2 = n1*n2 - U1;
U = min(U1,U2);
fprintf('p = %.4f, Mann-Whitney U = %.0f\n', p, U);


s = strcat("FigS GridPhaseDist_Group2.emf");
exportgraphics(fig2, s);


A = Y{1};
B = Y{2};
mA = mean(A);
mB = mean(B);
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');

fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
fprintf('Group2: %.3f ± %.3f\n', mB, sdB);





%% =================== Functions ===================

function Results = compute_within_between_phase(Gphase_table)

Results = [];

for m = 1:size(Gphase_table,1)
    for d = 1:size(Gphase_table,2)

        T = Gphase_table{m,d};
        if isempty(T)
            continue
        end

        T = T(T.CellID_1 ~= T.CellID_2,:);
        if isempty(T)
            continue
        end

        % module 1 & 2
        mod_ok = (T.ModuleInd_Cell1 <= 2) & (T.ModuleInd_Cell2 <= 2);

        % ===== within =====
        idx_w = (T.WithinOrBetween == 1) & mod_ok;
        phase_w = T.PhaseDistance(idx_w);
        phys_w  = T.PhysicalDistance(idx_w);

        % ===== between =====
        idx_b = (T.WithinOrBetween == 0) & mod_ok;
        phase_b = T.PhaseDistance(idx_b);
        phys_b  = T.PhysicalDistance(idx_b);

        if length(phase_w) < 3 || length(phase_b) < 3
            continue
        end

        out.within_phase = phase_w;
        out.within_phys  = phys_w;
        out.between_phase = phase_b;
        out.between_phys  = phys_b;

        out.within_mean = mean(phase_w,'omitnan');
        out.within_std  = std(phase_w,'omitnan');

        out.between_mean = mean(phase_b,'omitnan');
        out.between_std  = std(phase_b,'omitnan');

        Results = [Results; out];
    end
end
end


%%
function [p] = compare_groups(data1, data2)

data1 = data1(~isnan(data1));
data2 = data2(~isnan(data2));

if length(data1)>2 && length(data2)>2
    h1 = kstest((data1-mean(data1))/std(data1));
    h2 = kstest((data2-mean(data2))/std(data2));
else
    h1=1; h2=1;
end

if h1==0 && h2==0
    [~,p] = ttest2(data1,data2);
else
    p = ranksum(data1,data2);
end
end

function plot_bar_with_sem(data_cell, pval, labels, ylab)

nGroups = length(data_cell);
means = zeros(1,nGroups);
sems  = zeros(1,nGroups);

for i = 1:nGroups
    means(i) = mean(data_cell{i},'omitnan');
    sems(i)  = std(data_cell{i},'omitnan')/sqrt(sum(~isnan(data_cell{i})));
end

bar(means); hold on;
errorbar(1:nGroups, means, sems, 'k','LineStyle','none','LineWidth',1.5);

set(gca,'XTick',1:nGroups,'XTickLabel',labels)
ylabel(ylab)

ymax = max(means+sems)*1.1;
line([1 2],[ymax ymax],'Color','k','LineWidth',1.5)

if pval < 0.001
    text(1.5,ymax,'***','HorizontalAlignment','center')
elseif pval < 0.01
    text(1.5,ymax,'**','HorizontalAlignment','center')
elseif pval < 0.05
    text(1.5,ymax,'*','HorizontalAlignment','center')
else
    text(1.5,ymax,'n.s.','HorizontalAlignment','center')
end

box off
end


%% shuffle

% =================== Helper function ===================
function [mean_bin, sem_bin, all_bins] = compute_binned_phase_module(Results, edges, type)
nBin = length(edges)-1;
all_bins = cell(nBin,1);

for i = 1:length(Results)
    if strcmp(type,'within')
        x = Results(i).within_phys;
        y = Results(i).within_phase;
    elseif strcmp(type,'between')
        x = Results(i).between_phys;
        y = Results(i).between_phase;
    else
        error('type must be "within" or "between"');
    end

    [~,~,bin_idx] = histcounts(x, edges);
    [N,e,idx] = histcounts(x, edges);
    for b = 1:nBin
        all_bins{b} = [all_bins{b}; y(bin_idx==b)];
    end
end

mean_bin = nan(nBin,1);
sem_bin  = nan(nBin,1);

for b = 1:nBin
    data = all_bins{b};
    data = data(~isnan(data));
    if ~isempty(data)
        mean_bin(b) = mean(data);
        sem_bin(b)  = std(data)/sqrt(length(data));
    end
end
end


%%
function [mean_bin, sem_bin, all_bins] = compute_binned_phase_single(x, y, edges)
nBin = length(edges)-1;
all_bins = cell(nBin,1);

[~,~,bin_idx] = histcounts(x, edges);

for b = 1:nBin
    all_bins{b} = y(bin_idx==b);
end

mean_bin = nan(nBin,1);
sem_bin  = nan(nBin,1);

for b = 1:nBin
    data = all_bins{b};
    data = data(~isnan(data));
    if ~isempty(data)
        mean_bin(b) = mean(data);
        sem_bin(b)  = std(data)/sqrt(length(data));
    end
end
end


%%
function [p_Group1_vs_Group2_fdr, p_Group1_vs_shuf_fdr, p_Group2_vs_shuf_fdr] = ...
    func_fourttest(bin_centers, allbins_Group1_within, allbins_Group2_within, shuf_all)

p_Group1_vs_Group2  = nan(length(bin_centers),1);
p_Group1_vs_shuf  = nan(length(bin_centers),1);
p_Group2_vs_shuf = nan(length(bin_centers),1);

for b = 1:length(bin_centers)

    neg_real  = allbins_Group1_within{b};
    Group2_real = allbins_Group2_within{b};
    neg_shuf  = shuf_all{1}(b,:);
    Group2_shuf = shuf_all{2}(b,:);

    neg_real  = neg_real(~isnan(neg_real));
    Group2_real = Group2_real(~isnan(Group2_real));
    neg_shuf  = neg_shuf(~isnan(neg_shuf));
    Group2_shuf = Group2_shuf(~isnan(Group2_shuf));

    if isempty(neg_real) || isempty(Group2_real) || isempty(neg_shuf) || isempty(Group2_shuf)
        continue
    end

    [~, p_Group1_vs_Group2(b)]  = ttest2(neg_real, Group2_real, 'Vartype','unequal');
    [~, p_Group1_vs_shuf(b)]  = ttest2(neg_real, neg_shuf,  'Vartype','unequal');
    [~, p_Group2_vs_shuf(b)] = ttest2(Group2_real, Group2_shuf,'Vartype','unequal');
end

% ===== FDR =====
p_Group1_vs_Group2_fdr  = mafdr(p_Group1_vs_Group2,'BHFDR',true);
p_Group1_vs_shuf_fdr  = mafdr(p_Group1_vs_shuf,'BHFDR',true);
p_Group2_vs_shuf_fdr = mafdr(p_Group2_vs_shuf,'BHFDR',true);

end



%%
function c = func_phasedistance_bar(i, allbins_Group1_within, allbins_Group2_within, shuf_all)
figure; hold on;

wt        = allbins_Group1_within{i};
wt_shuf   = shuf_all{1}(i,:);
eb2       = allbins_Group2_within{i};
eb2_shuf  = shuf_all{2}(i,:);

data_cell = {wt, wt_shuf, eb2, eb2_shuf};

% ===== mean / sem =====
nGroups = length(data_cell);
means = zeros(1,nGroups);
sems  = zeros(1,nGroups);

for k = 1:nGroups
    d = data_cell{k};
    means(k) = mean(d,'omitnan');
    sems(k)  = std(d,'omitnan') / sqrt(sum(~isnan(d)));
end


xpos = [1 2   4 5];

b = bar(xpos, means, 0.8); hold on;
b.FaceColor = 'flat';

colors = [
    0 0.4 1;   % WT
    0.6 0.8 1;   % WT shuf）
    1.0 0.2 0.2;   % EB2
    1 0.6 0.6    % EB2 shuf
    ];

for k = 1:nGroups
    b.CData(k,:) = colors(k,:);
end

errorbar(xpos, means, sems, 'k', 'LineStyle','none','LineWidth',1.5);

set(gca,'XTick',[1.5 4.5], ...
    'XTickLabel',{'WT','EB2'})

ylabel('phase distance')

ymax = max(means+sems)*1.2;

box off
ylim([0 0.5])

%%
data = [wt; wt_shuf'; eb2; eb2_shuf'];
Group = [ ...
    repmat({'WT'}, length(wt),1);
    repmat({'WT_shuf'}, length(wt_shuf),1);
    repmat({'EB2'},    length(eb2),1);
    repmat({'EB2_shuf'},    length(eb2_shuf),1)];
[p, tbl, stats] = anova1(data, Group, 'off');
figure;
[c,~,~,gnames] = multcompare(stats, 'CType','tukey-kramer');

end