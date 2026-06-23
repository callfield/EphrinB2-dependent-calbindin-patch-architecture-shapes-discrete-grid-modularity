%% ========== Load data ==========
clear; close all hidden

GROUP1 = load("../Group1_CellTables.mat");
GROUP2 = load("../Group2_CellTables.mat");
GROUP1_T = GROUP1.GROUP1_T;
GROUP2_T = GROUP2.GROUP2_T;




%% Use only grid cells
GROUP1_T_modified = GROUP1_T;
GROUP2_T_modified = GROUP2_T;

for m = 1:numel(GROUP1_T_modified)
    if ~isempty(GROUP1_T_modified{m})
        idx3 = GROUP1_T_modified{m}.Grid_scale == 0;
        GROUP1_T_modified{m}(idx3, :) = [];
    end
end

for m = 1:numel(GROUP2_T_modified)
    if ~isempty(GROUP2_T_modified{m})
        idx3 = GROUP2_T_modified{m}.Grid_scale == 0;
        GROUP2_T_modified{m}(idx3, :) = [];
    end
end

GROUP1_T = GROUP1_T_modified;
GROUP2_T = GROUP2_T_modified;


N_grp1 = cellfun(@height , GROUP1_T);
N_grp1_all = sum(N_grp1, 'all');
N_grp2 = cellfun(@height , GROUP2_T);
N_grp2_all = sum(N_grp2, 'all');

%%

% Group1 (GROUP1_T)
neg_mean = nan(size(GROUP1_T));
for s = 1:size(GROUP1_T, 1)
    for t = 1:size(GROUP1_T, 2)
        T = GROUP1_T{s,t};
    
        if ~isempty(T)
            neg_mean(s,t) = mean(T.Grid_Score, 'omitnan');
        else
            neg_mean(s,t) = NaN;
        end
    end
end


% Group2 (GROUP2_T)
Group2_mean = nan(size(GROUP2_T));
for s = 1:size(GROUP2_T, 1)
    for t = 1:size(GROUP2_T, 2)
        T = GROUP2_T{s,t};
    
        if ~isempty(T)
            Group2_mean(s,t) = mean(T.Grid_Score, 'omitnan');
        else
            Group2_mean(s,t) = NaN;
        end
    end
end

neg_mean = mean(neg_mean, 2);
Group2_mean = mean(Group2_mean, 2);


%%
% CGroup1 = [1 1 1];
% CGroup2 = [255 168 60]/255;
CGroup1 = [1 1 1];
CGroup2 = [240 134 134]/255;

pos = [300 500 130 160];
fig = figure('Position', pos);

neg = neg_mean;
Group2 = Group2_mean;
means = [
    mean(neg,'omitnan')
    mean(Group2,'omitnan')
    ];

sems = [
    std(neg,'omitnan')/sqrt(sum(~isnan(neg)))
    std(Group2,'omitnan')/sqrt(sum(~isnan(Group2)))
    ];

x = [0.9, 2.1];
b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
hold on
b.CData(1,:) = [1 1 1];
b.CData(2,:) = CGroup2;

errorbar(x, means, sems, ...
    'k','LineStyle','none','LineWidth',0.5,'CapSize',12)

ytol = 0.1;
dx = 0.15;
x1 = SimpleBeeSwarm(neg, x(1), ytol, dx);
x2 = SimpleBeeSwarm(Group2, x(2), ytol, dx);
scatter(x1, neg, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
scatter(x2, Group2, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)

xlim([0 3])
ylim([0 4])
yticks([0 2 4])

set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{})

box off


A = neg(~isnan(neg));
B = Group2(~isnan(Group2));

[~, p_ttest, ~, stats] = ttest2(A, B, 'Vartype','unequal');
tval = stats.tstat;
df   = stats.df;
fprintf('p = %.4f, t(%.2f) = %.3f\n', p_ttest, df, tval);

s = strcat("GridScore.emf");
exportgraphics(fig, s);



%% functions
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