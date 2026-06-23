clear all;

addpath('function')

load("../Data.mat")

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

%% Grid scale scatter with KSD regression slope
gp1 = [];
for s = 1:7
    for t = 1:3
        data = GROUP1{s,t}(GROUP1{s,t}(:,3) > 0, 1:13);
        gp1 = [gp1; data];
    end
end

gp2 = [];
for s = 1:7
    for t = 1:3
        data = GROUP2{s,t}(GROUP2{s,t}(:,3) > 0, 1:13);
        gp2 = [gp2; data];
    end
end

data1 = gp1;
data2 = gp2;
peaks1 = ksdnorm200(data1);
peaks2 = ksdnorm200(data2);

close all
figure('Position', [100, 100, 250, 500]);
xlims = [25 90];
ylims = [0 800];

plot_GridScaleVsLocation_KSDRegression_2(data1, peaks1, xlims, ylims);
exportgraphics(gcf, strcat("GScale_ScatterwithSlope gp1 all trial.pdf"))
close all

figure('Position', [100, 100, 250, 500]);
plot_GridScaleVsLocation_KSDRegression_2(data2, peaks2, xlims, ylims);
exportgraphics(gcf, strcat("GScale_ScatterwithSlope gp2 all trial.pdf"))
close all

%% Region-adjusted Grid Scale distribution
% data1: group 1
% data2: group 2
% Column 1: sampling region, e.g. DV position
% Column 3: Grid Scale

GS_COL = 3;
REGION_COL = 1;

group1Name = 'group1';
group2Name = 'group2';

bandwid = 3.5;
nRegionBins = 8;
nPerm = 5000;

xRange = [20 100];
x = linspace(xRange(1), xRange(2), 300);

binningMode = "quantile";
% "quantile"  : make sampling-region bins with approximately equal sample counts.
% "equalwidth": make sampling-region bins with equal physical widths.

rng(1);  % Reproducible permutation test.

%% Extract variables
gs1 = data1(:, GS_COL);
gs2 = data2(:, GS_COL);

reg1 = data1(:, REGION_COL);
reg2 = data2(:, REGION_COL);

% Remove NaN and Inf values.
valid1 = isfinite(gs1) & isfinite(reg1);
valid2 = isfinite(gs2) & isfinite(reg2);

gs1 = gs1(valid1);
gs2 = gs2(valid2);
reg1 = reg1(valid1);
reg2 = reg2(valid2);

%% Bin sampling region
allReg = [reg1; reg2];

switch binningMode
    case "quantile"
        edges = quantile(allReg, linspace(0, 1, nRegionBins + 1));
        edges = unique(edges, 'stable');

    case "equalwidth"
        edges = linspace(min(allReg), max(allReg), nRegionBins + 1);

    otherwise
        error('binningMode must be "quantile" or "equalwidth".');
end

% If the sampling region is nearly constant, use one catch-all bin.
if numel(edges) < 2
    edges = [-inf inf];
else
    edges(1) = -inf;
    edges(end) = inf;
end

bin1 = discretize(reg1, edges);
bin2 = discretize(reg2, edges);

ok1 = ~isnan(bin1);
ok2 = ~isnan(bin2);

gs1 = gs1(ok1);
gs2 = gs2(ok2);
reg1 = reg1(ok1);
reg2 = reg2(ok2);
bin1 = bin1(ok1);
bin2 = bin2(ok2);

nBin = numel(edges) - 1;

%% Count samples per sampling-region bin
c1 = accumarray(bin1, 1, [nBin 1], @sum, 0);
c2 = accumarray(bin2, 1, [nBin 1], @sum, 0);

commonBin = c1 > 0 & c2 > 0;

if ~any(commonBin)
    error('No overlapping sampling-region bins between the two groups.');
end

% Use only sampling-region bins common to both groups.
% The target distribution is the per-bin minimum count, excluding bins present in only one group.
target = min(c1, c2);
target(~commonBin) = 0;
target = target / sum(target);

%% Compute region-adjusted weights
w1 = zeros(size(gs1));
w2 = zeros(size(gs2));

for b = 1:nBin
    if target(b) > 0
        w1(bin1 == b) = target(b) / c1(b);
        w2(bin2 == b) = target(b) / c2(b);
    end
end

% Remove data outside the common sampling-region support.
keep1 = w1 > 0;
keep2 = w2 > 0;

gs1 = gs1(keep1);
gs2 = gs2(keep2);
reg1 = reg1(keep1);
reg2 = reg2(keep2);
bin1 = bin1(keep1);
bin2 = bin2(keep2);
w1 = w1(keep1);
w2 = w2(keep2);

% Normalize weights.
w1 = w1 / sum(w1);
w2 = w2 / sum(w2);

%% Region-adjusted KDE
f1 = ksdensity(gs1, x, ...
    'Bandwidth', bandwid, ...
    'Weights', w1);

f2 = ksdensity(gs2, x, ...
    'Bandwidth', bandwid, ...
    'Weights', w2);

%% Region-adjusted KS-like statistic
D_adj = weightedKS(gs1, gs2, w1, w2);

%% Region-stratified permutation test
% Shuffle group labels only within each sampling-region bin.
allGS = [gs1; gs2];
allBin = [bin1; bin2];

label = [ones(size(gs1)); 2 * ones(size(gs2))];

D_perm = nan(nPerm, 1);
validBins = find(target > 0);

for ip = 1:nPerm

    permLabel = label;

    for b = validBins'
        idx = find(allBin == b);

        % Preserve the number of group 1 samples within each sampling-region bin.
        n1b = sum(label(idx) == 1);

        idxShuffle = idx(randperm(numel(idx)));

        permLabel(idx) = 2;
        permLabel(idxShuffle(1:n1b)) = 1;
    end

    % Recompute weights after permutation.
    wp = zeros(size(allGS));

    for b = validBins'
        idx1 = allBin == b & permLabel == 1;
        idx2 = allBin == b & permLabel == 2;

        wp(idx1) = target(b) / sum(idx1);
        wp(idx2) = target(b) / sum(idx2);
    end

    pgs1 = allGS(permLabel == 1);
    pgs2 = allGS(permLabel == 2);

    pw1 = wp(permLabel == 1);
    pw2 = wp(permLabel == 2);

    pw1 = pw1 / sum(pw1);
    pw2 = pw2 / sum(pw2);

    D_perm(ip) = weightedKS(pgs1, pgs2, pw1, pw2);
end

p_perm = (sum(D_perm >= D_adj) + 1) / (nPerm + 1);
H_adj = p_perm < 0.05;

%% Plot
fig = figure('Position', [100, 100, 500, 200]);
hold on

h1 = area(x, f1, ...
    'FaceColor', [0.3 0.3 0.3], ...
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');

h2 = area(x, f2, ...
    'FaceColor', [0.9 0.3 0.3], ...
    'FaceAlpha', 0.5, ...
    'EdgeColor', 'none');

plot(x, f1, ...
    'Color', [0.3 0.3 0.3] * 0.8, ...
    'LineWidth', 2);

plot(x, f2, ...
    'Color', [0.9 0.3 0.3] * 0.8, ...
    'LineWidth', 2);

xlim(xRange)
yticks([])

xlabel('Grid Scale')
ylabel('Region-adjusted density')

legend([h1 h2], {group1Name, group2Name}, ...
    'Location', 'best', ...
    'Box', 'off');

title(sprintf('Region-adjusted KS-like D = %.3f, p_{perm} = %.4f', ...
    D_adj, p_perm), ...
    'Interpreter', 'tex');

set(gca, 'FontSize', 11);

%% Export
outFile = "G_scale_regionAdjusted_all_trial.emf";
exportgraphics(fig, outFile, "ContentType", "vector");

close

[H, p, ksstat] = kstest2(gs1, gs2)

function D = weightedKS(x1, x2, w1, w2)
% KS-like statistic based on weighted empirical CDFs.

    x1 = x1(:);
    x2 = x2(:);
    w1 = w1(:);
    w2 = w2(:);

    w1 = w1 / sum(w1);
    w2 = w2 / sum(w2);

    allx = unique(sort([x1; x2]));

    [u1, ~, ic1] = unique(x1);
    mass1 = accumarray(ic1, w1, [numel(u1), 1], @sum, 0);
    cdf1 = cumsum(mass1);

    [u2, ~, ic2] = unique(x2);
    mass2 = accumarray(ic2, w2, [numel(u2), 1], @sum, 0);
    cdf2 = cumsum(mass2);

    F1 = stepCDF(u1, cdf1, allx);
    F2 = stepCDF(u2, cdf2, allx);

    D = max(abs(F1 - F2));
end

function F = stepCDF(u, cdfval, x)
% Helper function to evaluate a stepwise CDF.

    u = u(:);
    cdfval = cdfval(:);
    x = x(:);

    if isscalar(u)
        F = double(x >= u(1)) * cdfval(end);
    else
        F = interp1(u, cdfval, x, 'previous', 'extrap');
        F(x < u(1)) = 0;
        F(x >= u(end)) = cdfval(end);
    end
end
