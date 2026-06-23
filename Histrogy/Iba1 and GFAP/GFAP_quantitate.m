%% ========== Load data ==========
clear; close all hidden

% read csv  	1 slice, 2 Area,	3 intensity(Mean), 4 intensity(Max)			5, Area

% negative control
N{1} = readmatrix('Nega1.csv');
N{2} = readmatrix('Nega2.csv');
N{3} = readmatrix('Nega3.csv');

% Casp(Cre positive)
P{1} = readmatrix('Casp1.csv');
P{2} = readmatrix('Casp2.csv');
P{3} = readmatrix('Casp3.csv');
P{4} = readmatrix('Casp4.csv');


%%

NormalI_N = cell(1,3);
for i = 1:3
    temp = N{i};
    temp(1,:) = [];
    NormalI_N{i} = nan(size(temp, 1)/2, 1);
    for j = 1:size(temp, 1)/2
        intensity_temp = temp(1 + 2*(j-1), 3);
        reference_intensity = temp(2*j, 3);
        NormalI_N{i}(j,1) = intensity_temp/reference_intensity;
    end
end

NormalI_P = cell(1,4);
for i = 1:4
    temp = P{i};
    temp(1,:) = [];
    NormalI_P{i} = nan(size(temp, 1)/2, 1);
    for j = 1:size(temp, 1)/2
        intensity_temp = temp(1 + 2*(j-1), 3);
        reference_intensity = temp(2*j, 3);
        NormalI_P{i}(j,1) = intensity_temp/reference_intensity;
    end
end

%%
NormalI_N_mean = cell2mat(cellfun(@(x) mean(x), NormalI_N, 'UniformOutput',false));
NormalI_P_mean = cell2mat(cellfun(@(x) mean(x), NormalI_P, 'UniformOutput',false));

%% animal-wise

Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

pos = [300 500 130 160];
fig = figure('Position', pos);


neg = NormalI_N_mean / mean(NormalI_N_mean,'omitnan');
casp = NormalI_P_mean / mean(NormalI_N_mean,'omitnan');


means = [
    mean(neg,'omitnan')
    mean(casp,'omitnan')
    ];
sems = [
    std(neg,'omitnan')/sqrt(sum(~isnan(neg)))
    std(casp,'omitnan')/sqrt(sum(~isnan(casp)))
    ];

x = [0.9, 2.1];
b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
hold on

b.CData(1,:) = [1 1 1];
b.CData(2,:) = Ccasp;

errorbar(x, means, sems, ...
    'k','LineStyle','none','LineWidth',0.5,'CapSize',12)


ytol = 0.05;
dx = 0.15;
x1 = SimpleBeeSwarm(neg, x(1), ytol, dx);
x2 = SimpleBeeSwarm(casp, x(2), ytol, dx);
scatter(x1, neg, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
scatter(x2, casp, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)


xlim([0 3])
ylim([0 2])
yticks([0 1 2])

set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{})

box off

A = neg(~isnan(neg));
B = casp(~isnan(casp));

[~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','equal');
tval = stats.tstat;
df   = stats.df;
fprintf('p = %.4f, t(%d) = %.3f\n', p_ttest, df, tval);


% [~, p_ttest, ~, stats] = ttest2(A, B, 'Vartype','unequal');
% tval = stats.tstat;
% df   = stats.df;
% fprintf('p = %.4f, t(%.2f) = %.3f\n', p_ttest, df, tval);

% s = strcat("GFAP_normalIntensity_Animal.emf");
% exportgraphics(fig, s);


%% slice-wise

Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

pos = [300 500 130 160];
fig = figure('Position', pos);


neg = cell2mat(NormalI_N(:)) / mean(cell2mat(NormalI_N(:)),'omitnan');
casp = cell2mat(NormalI_P(:)) / mean(cell2mat(NormalI_N(:)),'omitnan');

means = [
    mean(neg,'omitnan')
    mean(casp,'omitnan')
    ];

sems = [
    std(neg,'omitnan')/sqrt(sum(~isnan(neg)))
    std(casp,'omitnan')/sqrt(sum(~isnan(casp)))
    ];

x = [0.9, 2.1];
b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
hold on

b.CData(1,:) = [1 1 1];
b.CData(2,:) = Ccasp;

errorbar(x, means, sems, ...
    'k','LineStyle','none','LineWidth',0.5,'CapSize',12)


ytol = 0.05;
dx = 0.15;
x1 = SimpleBeeSwarm(neg, x(1), ytol, dx);
x2 = SimpleBeeSwarm(casp, x(2), ytol, dx);
scatter(x1, neg, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
scatter(x2, casp, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)


xlim([0 3])
ylim([0 2.5])
yticks([0 1 2 3])

set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{})

box off

A = neg(~isnan(neg));
B = casp(~isnan(casp));

% [~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','equal');
[~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','unequal');
tval = stats.tstat;
df   = stats.df;
fprintf('p = %.4f, t(%d) = %.3f\n', p_ttest, df, tval);


s = strcat("GFAP_normalIntensity_Slice.emf");
exportgraphics(fig, s);

mA = mean(A);
mB = mean(B);
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');

fprintf('Control: %.3f ± %.3f\n', mA, sdA);
fprintf('Casp: %.3f ± %.3f\n', mB, sdB);


%% functions
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

