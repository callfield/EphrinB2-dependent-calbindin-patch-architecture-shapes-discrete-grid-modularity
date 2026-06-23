%% ========== Load data ==========
clear; close all hidden


GROUP1 = load("../Group1_CellTables.mat");
GROUP2 = load("../Group2_CellTables.mat");
GROUP1_T = GROUP1.GROUP1_T;
GROUP2_T = GROUP2.GROUP2_T;



%% load spike timing
pop_str = {'GridCell', 'NonGridSpatialCell', 'AllCell'};
group_animal_str = {'Ctrl', 'Exp'};


for pop = 1 %1:3 
    fprintf('\n=============================\n');
    fprintf('POPULATION: %d / 3\n', pop);
    fprintf('=============================\n');

    FindGrid = pop;

    [GROUP1_T, GROUP2_T, keepIdxgrp1a, keepIdxgrp2] = ...
        func_findGridOrNonGrid(GROUP1_T, GROUP2_T, FindGrid);



    %%
    Frate_results = cell(2,1);

    for group_animal = 1:2
        if group_animal == 1
            OcT_temp = GROUP1_T;
            group_name = 'Group1';
        else
            OcT_temp = GROUP2_T;
            group_name = 'Group2';
        end
        fprintf('\n  >> Group: %s (%d/2)\n', group_name, group_animal);

        Max_m = size(OcT_temp, 1);
        Max_d = size(OcT_temp, 2);
        Frate_results{group_animal} = cell(Max_m, Max_d);
        for m = 1:Max_m
            for d = 1:Max_d
                fprintf('    Mouse %d/%d | Day %d/%d\n', m, Max_m, d, Max_d);
                close all

                OcT = OcT_temp{m, d};
                if isempty(OcT)
                    continue;
                end
                
                    % calcium traces
                    A = OcT.dF_source;
                    filePaths = cellfun(@(x) x.Properties.Source, A, 'UniformOutput',false);
       
                    new_root = '..\01_Each_animal';

                    filePaths{i} = strrep(filePaths{i}, new_root);


                    % loaded_CaTrace = load(filePaths{1});
                    [filepath,~,~] = fileparts(filePaths{1});

                    % load animal tracks and speeds
                    A = OcT.Trk_source;
                    filePaths_track = cellfun(@(x) x.Properties.Source, A, 'UniformOutput',false);
                    [~,name,ext] = fileparts(filePaths_track{1});
                    fname_track = fullfile(filepath, strcat(name, ext));
                    loaded_AnimalTrack = load(fname_track);

                    AnimalTrack = loaded_AnimalTrack.Trk_withTimeStamp; % timestamp already aligned
                    AnimalTrack(:,1) = AnimalTrack(:,1) - AnimalTrack(1,1); % set start = zero

                    % load spike timing
                    A = load(fullfile(filepath, "ST_dF_grid_aut_data.mat"));
                    SpkInd = A.m_lk;

                    if group_animal == 1
                        ind = keepIdxgrp1a{m, d}; 
                    elseif group_animal == 2
                        ind = keepIdxgrp2{m, d};
                    end
                    SpkInd_selectedCells = SpkInd(ind);
                    dt = 0.1; % sampling interval
                    Frate = cellfun(@(x) length(x)/(size(AnimalTrack,1)*dt), ...
                        SpkInd_selectedCells);

                    Frate_results{group_animal}{m,d} = Frate;
            end
        end
    end


    %% ===== session-wise mean =====

    % CGroup1 = [1 1 1];
    % Cgroup2 = [255 168 60]/255;

    CGroup1 = [1 1 1];
    Cgroup2 = [240 134 134]/255;


   Group1_Frate = Frate_results{1};
   Group1_Frate = cellfun(@(x) mean(x, 'omitmissing'),Group1_Frate);

    Group2_Frate = Frate_results{2};
    Group2_Frate = cellfun(@(x) mean(x, 'omitmissing'), Group2_Frate);


    pos = [300 500 130 160];
    fig = figure('Position', pos);


    grp1 =Group1_Frate(:);
    grp2 = Group2_Frate(:);

    means = [
        mean(grp1,'omitnan')
        mean(grp2,'omitnan')
        ];

    sems = [
        std(grp1,'omitnan')/sqrt(sum(~isnan(grp1)))
        std(grp2,'omitnan')/sqrt(sum(~isnan(grp2)))
        ];

    x = [0.9, 2.1];
    b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
    hold on

    b.CData(1,:) = [1 1 1];
    b.CData(2,:) = Cgroup2;

    errorbar(x, means, sems, ...
        'k','LineStyle','none','LineWidth',0.5,'CapSize',12)


    ytol = 0.002;
    dx = 0.15;
    x1 = SimpleBeeSwarm(grp1, x(1), ytol, dx);
    x2 = SimpleBeeSwarm(grp2, x(2), ytol, dx);
    scatter(x1, grp1, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
    scatter(x2, grp2, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)

    xlim([0 3])
    ylim([0 0.2])
    yticks([0 0.1 0.2])

    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{})

    box off

    A = grp1(~isnan(grp1));
    B = grp2(~isnan(grp2));

    [~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','equal');
    tval = stats.tstat;
    df   = stats.df;
    fprintf('%s  p = %.4f, t(%d) = %.3f\n', pop_str{pop} , p_ttest, df, tval);


    % [~, p_ttest, ~, stats] = ttest2(A, B, 'Vartype','unequal');
    % tval = stats.tstat;
    % df   = stats.df;
    % fprintf('p = %.4f, t(%.2f) = %.3f\n', p_ttest, df, tval);


    % s = strcat("SessionWise_Frate", pop_str{pop}, '.emf');
    % exportgraphics(fig, s);

        %% ===== Animal-wise analysis =====

   Group1_animal = cell(size(GROUP1_T,1),1);
    Group2_animal = cell(size(GROUP2_T,1),1);

    % Ctrl
    for m = 1:size(Frate_results{1},1)
        tmp = Frate_results{1}(m,:);
        tmp = cellfun(@(x) mean(x,'omitnan'), tmp);
       Group1_animal{m} = mean(tmp,'omitnan');
    end

    % grp2
    for m = 1:size(Frate_results{2},1)
        tmp = Frate_results{2}(m,:);
        tmp = cellfun(@(x) mean(x,'omitnan'), tmp);
        Group2_animal{m} = mean(tmp,'omitnan');
    end

    grp1 = cell2mat(Ctrl_animal);
    grp2 = cell2mat(Group2_animal);

    
    pos = [300 500 130 160];
    fig = figure('Position', pos);

    means = [
        mean(grp1,'omitnan')
        mean(grp2,'omitnan')
        ];

    sems = [
        std(grp1,'omitnan')/sqrt(sum(~isnan(grp1)))
        std(grp2,'omitnan')/sqrt(sum(~isnan(grp2)))
        ];

    x = [0.9, 2.1];
    b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
    hold on

    b.CData(1,:) = [1 1 1];
    b.CData(2,:) = Cgroup2;

    errorbar(x, means, sems, ...
        'k','LineStyle','none','LineWidth',0.5,'CapSize',12)


    ytol = 0.002;
    dx = 0.15;
    x1 = SimpleBeeSwarm(grp1, x(1), ytol, dx);
    x2 = SimpleBeeSwarm(grp2, x(2), ytol, dx);
    scatter(x1, grp1, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
    scatter(x2, grp2, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)

    % xlim([0.5 2.5])
    xlim([0 3])
    ylim([0 0.2])
    yticks([0 0.1 0.2])

    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{})

    box off

    A = grp1(~isnan(grp1));
    B = grp2(~isnan(grp2));

    % [~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','equal');
    [~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','unequal');
    tval = stats.tstat;
    df   = stats.df;
    fprintf('%s  p = %.4f, t(%d) = %.3f\n', pop_str{pop} , p_ttest, df, tval);


    s = strcat("AnimalWise_Frate ", pop_str{pop}, '.emf');
    exportgraphics(fig, s);


    mA = mean(A);
mB = mean(B);
sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Group1: %.3f ± %.3f\n', mA, sdA);
fprintf('Group2 %.3f ± %.3f\n', mB, sdB);

    %% Cell-Wise analysis


   Group1_Frate = cell2mat(Frate_results{1}(:));
    Group2_Frate =  cell2mat(Frate_results{2}(:));

    pos = [500 500 130 160];
    fig = figure('Position', pos);

    grp1 =Group1_Frate(:);
    grp2 = Group2_Frate(:);

    means = [
        mean(grp1,'omitnan')
        mean(grp2,'omitnan')
        ];

    sems = [
        std(grp1,'omitnan')/sqrt(sum(~isnan(grp1)))
        std(grp2,'omitnan')/sqrt(sum(~isnan(grp2)))
        ];

    x = [0.9, 2.1];
    b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
    hold on

    b.CData(1,:) = [1 1 1];
    b.CData(2,:) = Cgroup2;

    errorbar(x, means, sems, ...
        'k','LineStyle','none','LineWidth',0.5,'CapSize',12)


    xlim([0 3])
    ylim([0 0.2])
    %% 
    yticks([0 0.1 0.2])

    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{})

    box off

    A = grp1(~isnan(grp1));
    B = grp2(~isnan(grp2));

    [~, p_ttest, ~, stats] = ttest2(A,B, 'Vartype','equal');
    tval = stats.tstat;
    df   = stats.df;
    fprintf('%s  p = %.4f, t(%d) = %.3f\n', pop_str{pop} , p_ttest, df, tval);

    % [~, p_ttest, ~, stats] = ttest2(A, B, 'Vartype','unequal');
    % tval = stats.tstat;
    % df   = stats.df;
    % fprintf('p = %.4f, t(%.2f) = %.3f\n', p_ttest, df, tval);

    % s = strcat("CellWise_Frate", pop_str{pop}, '.emf');
    % exportgraphics(fig, s);


end





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


%% grid cell detection

function [GROUP1_T, GROUP2_T, keepIdxgrp1a, keepIdxgrp2] = func_findGridOrNonGrid(GROUP1_T, GROUP2_T, FindGrid)

GROUP1_T = GROUP1_T;
GROUP2_T = GROUP2_T;

keepIdxgrp1a = cell(size(GROUP1_T));
keepIdxgrp2 = cell(size(GROUP2_T));

for m = 1:numel(GROUP1_T)

    if isempty(GROUP1_T{m})
        continue
    end

    switch FindGrid
        case 1 % grid
            keepIdx = GROUP1_T{m}.Grid_scale > 0;

        case 2 % non-grid spatial
            keepIdx = GROUP1_T{m}.Grid_scale == 0 & ...
                      GROUP1_T{m}.Spatial_Score >= 10;

        case 3 % all
            keepIdx = true(height(GROUP1_T{m}),1);
    end

    keepIdxgrp1a{m} = find(keepIdx);
    GROUP1_T{m} = GROUP1_T{m}(keepIdx,:);
end


for m = 1:numel(GROUP2_T)

    if isempty(GROUP2_T{m})
        continue
    end

    switch FindGrid
        case 1 % grid
            keepIdx = GROUP2_T{m}.Grid_scale > 0;

        case 2 % non-grid spatial
            keepIdx = GROUP2_T{m}.Grid_scale == 0 & ...
                      GROUP2_T{m}.Spatial_Score >= 10;

        case 3 % all
            keepIdx = true(height(GROUP2_T{m}),1);
    end

    keepIdxgrp2{m} = find(keepIdx);
    GROUP2_T{m} = GROUP2_T{m}(keepIdx,:);
end

end
