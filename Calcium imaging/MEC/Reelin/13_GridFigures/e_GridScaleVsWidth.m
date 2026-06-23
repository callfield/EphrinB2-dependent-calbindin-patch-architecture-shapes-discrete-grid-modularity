%% ========== Load data ==========
clear; close all hidden


GROUP1 = load("../Group1_CellTables.mat");
GROUP2 = load("../Group2_CellTables.mat");
GROUP1_T = GROUP1.GROUP1_T;
GROUP2_T = GROUP2.GROUP2_T;



%% ==================

% CGROUP1 = [1 1 1];
% CGROUP2 = [255 168 60]/255;
CGROUP1 = [1 1 1];
CGROUP2 = [240 134 134]/255;

% close all


%% =========================================================
% Grid scale vs Grid field width
%% =========================================================

for TrialSelection = 0 %[0,1]

    close all


    %% =========================================================
    % GROUP1
    %% =========================================================

    R_gp1 = [];
    Slope_gp1 = [];

    for m = 1:size(GROUP1,1)
        for d = 1:size(GROUP1,2)

            T = GROUP1{m,d};

            if isempty(T) || height(T)==0
                continue
            end

            if TrialSelection == 1
                if MaxModule_gp1(m,d) == 0
                    continue
                end
            end

            % ===== use all cells =====
            idx = ...
                (T.Grid_scale > 0) & ...
                (T.Grid_field > 0) & ...
                (~isnan(T.Grid_scale)) & ...
                (~isnan(T.Grid_field));

            x = T.Grid_scale(idx);
            y = T.Grid_field(idx);

            if length(x) < 5
                continue
            end

            % ===== correlation =====
            [R,P] = corr(x, y, 'rows','complete');

            % ===== linear fit =====
            p = polyfit(x, y, 1);
            slope = p(1);

            R_gp1(end+1,1) = R;
            Slope_gp1(end+1,1) = slope;

        end
    end


    %% =========================================================
    % GROUP2
    %% =========================================================

    R_gp2 = [];
    Slope_gp2 = [];

    for m = 1:size(GROUP2_T,1)
        for d = 1:size(GROUP2_T,2)

            T = GROUP2_T{m,d};

            if isempty(T) || height(T)==0
                continue
            end

            if TrialSelection == 1
                if MaxModule_Exp(m,d) == 0
                    continue
                end
            end

            idx = ...
                (T.Grid_scale > 0) & ...
                (T.Grid_field > 0) & ...
                (~isnan(T.Grid_scale)) & ...
                (~isnan(T.Grid_field));

            x = T.Grid_scale(idx);
            y = T.Grid_field(idx);

            if length(x) < 5
                continue
            end

            % ===== correlation =====
            [R,P] = corr(x, y, 'rows','complete');

            % ===== linear fit =====
            p = polyfit(x, y, 1);
            slope = p(1);

            R_gp2(end+1,1) = R;
            Slope_gp2(end+1,1) = slope;

        end
    end


    %% =========================================================
    % statistics
    %% =========================================================

    [~, p_R] = ttest2(R_gp1, R_gp2);
    [~, p_slope] = ttest2(Slope_gp1, Slope_gp2);

    fprintf('\n===== Correlation R =====\n')
    fprintf('GROUP1 : %.3f ± %.3f\n', ...
        mean(R_gp1,'omitnan'), ...
        std(R_gp1,'omitnan')/sqrt(sum(~isnan(R_gp1))))

    fprintf('GROUP2    : %.3f ± %.3f\n', ...
        mean(R_gp2,'omitnan'), ...
        std(R_gp2,'omitnan')/sqrt(sum(~isnan(R_gp2))))

    fprintf('p = %.4f\n', p_R)


    fprintf('\n===== Slope =====\n')
    fprintf('GROUP1 : %.3f ± %.3f\n', ...
        mean(Slope_gp1,'omitnan'), ...
        std(Slope_gp1,'omitnan')/sqrt(sum(~isnan(Slope_gp1))))

    fprintf('GROUP2    : %.3f ± %.3f\n', ...
        mean(Slope_gp2,'omitnan'), ...
        std(Slope_gp2,'omitnan')/sqrt(sum(~isnan(Slope_gp2))))

    fprintf('p = %.4f\n', p_slope)


    %% =========================================================
    % Figure : Correlation R
    %% =========================================================

    % fig = figure('Position',[300 500 120 140]);
    pos = [300 500 130 160];
    fig = figure('Position', pos);

    means = [
        mean(R_gp1,'omitnan')
        mean(R_gp2,'omitnan')
        ];

    sems = [
        std(R_gp1,'omitnan')/sqrt(sum(~isnan(R_gp1)))
        std(R_gp2,'omitnan')/sqrt(sum(~isnan(R_gp2)))
        ];

    x = [0.9, 2.1];
    b = bar(x, means,'FaceColor','flat', 'BarWidth',0.6);
    hold on

    b.CData(1,:) = [1 1 1];
    b.CData(2,:) = CGROUP2;

    errorbar(x, means, sems, ...
        'k','LineStyle','none','LineWidth',0.5)


    ytol = 0.03;
    dx = 0.15;
    x1 = SimpleBeeSwarm(R_gp1, x(1), ytol, dx);
    x2 = SimpleBeeSwarm(R_gp2, x(2), ytol, dx);
    scatter(x1, R_gp1, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)
    scatter(x2, R_gp2, 8, 'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth',0.5)

    % xlim([0.5 2.5])
    xlim([0 3])

    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{})

    ylabel('R')

    box off



    yticks([-0.2 0 0.5 1.0])
    ylim([-0.2 1])
    ylabel('')
    xticks([])
    yticklabels([])




    %% =========================================================
    % pooled data
    %% =========================================================

    AllScale_gp1 = [];
    AllField_gp1 = [];

    AllScale_gp2 = [];
    AllField_gp2 = [];

    %% =========================================================
    % EACH SESSION : GROUP1
    %% =========================================================

    fig = figure('Position',[100 100 900 500]);

    plot_i = 1;

    for m = 1:size(GROUP1,1)
        for d = 1:size(GROUP1,2)

            T = GROUP1{m,d};

            if isempty(T) || height(T)==0
                continue
            end

            if TrialSelection == 1
                if MaxModule_gp1(m,d) == 0
                    continue
                end
            end

            idx = ...
                (T.Grid_scale > 0) & ...
                (T.Grid_field > 0) & ...
                (~isnan(T.Grid_scale)) & ...
                (~isnan(T.Grid_field));

            x = T.Grid_scale(idx);
            y = T.Grid_field(idx);

            % if length(x) < 5
            %     continue
            % end

            % pooled data
            AllScale_gp1 = [AllScale_gp1; x];
            AllField_gp1 = [AllField_gp1; y];

            % ===== subplot =====
            subplot(3,5,plot_i)
            hold on

            scatter(x, y, 10, ...
                'MarkerFaceColor','w', ...
                'MarkerEdgeColor',CGROUP1, ...
                'LineWidth',0.5)

            % linear fit
            p = polyfit(x,y,1);

            xx = linspace(min(x), max(x), 100);
            yy = polyval(p, xx);

            plot(xx, yy, 'k', 'LineWidth',1)

            % correlation
            [R,P] = corr(x,y,'rows','complete');

            title(sprintf('M%d D%d\nR=%.2f',m,d,R), ...
                'FontSize',8)

            box off

            xlabel('Scale')
            ylabel('Field')

            plot_i = plot_i + 1;

        end
    end

    sgtitle('GROUP1 : each session')


    %% =========================================================
    % EACH SESSION : GROUP2
    %% =========================================================

    fig = figure('Position',[100 100 900 500]);

    plot_i = 1;

    for m = 1:size(GROUP2_T,1)
        for d = 1:size(GROUP2_T,2)

            T = GROUP2_T{m,d};

            if isempty(T) || height(T)==0
                continue
            end

            if TrialSelection == 1
                if MaxModule_Exp(m,d) == 0
                    continue
                end
            end

            idx = ...
                (T.Grid_scale > 0) & ...
                (T.Grid_field > 0) & ...
                (~isnan(T.Grid_scale)) & ...
                (~isnan(T.Grid_field));

            x = T.Grid_scale(idx);
            y = T.Grid_field(idx);

            % if length(x) < 5
            %     continue
            % end

            % pooled data
            AllScale_gp2 = [AllScale_gp2; x];
            AllField_gp2 = [AllField_gp2; y];

            % ===== subplot =====
            subplot(3,7,plot_i)
            hold on

            scatter(x, y, 10, ...
                'MarkerFaceColor','w', ...
                'MarkerEdgeColor',CGROUP2, ...
                'LineWidth',0.5)

            % linear fit
            p = polyfit(x,y,1);

            xx = linspace(min(x), max(x), 100);
            yy = polyval(p, xx);

            plot(xx, yy, 'Color',CGROUP2*0.8, 'LineWidth',1)

            % correlation
            [R,P] = corr(x,y,'rows','complete');

            title(sprintf('M%d D%d\nR=%.2f',m,d,R), ...
                'FontSize',8)

            box off

            xlabel('Scale')
            ylabel('Field')

            plot_i = plot_i + 1;

        end
    end

    sgtitle('GROUP2 : each session')


    %% =========================================================
    % POOLED ALL SESSIONS
    %% =========================================================

    % fig = figure('Position',[300 300 250 220]);
    % hold on

    pos2 = [300   500   160   160];
    fig = figure('Position',pos2);
    hold on

    % ===== GROUP1 =====
    scatter(AllScale_gp1, AllField_gp1, ...
        10, ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.5)

    % fit
    p1 = polyfit(AllScale_gp1, AllField_gp1, 1);

    xx1 = linspace(min(AllScale_gp1), max(AllScale_gp1), 100);
    yy1 = polyval(p1, xx1);

    plot(xx1, yy1, 'k', 'LineWidth',2)

    % correlation
    [R1,P1] = corr(AllScale_gp1, AllField_gp1, ...
        'rows','complete');


    % ===== GROUP2 =====
    scatter(AllScale_gp2, AllField_gp2, ...
        10, ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor',CGROUP2, ...
        'LineWidth',0.5)

    % fit
    p2 = polyfit(AllScale_gp2, AllField_gp2, 1);

    xx2 = linspace(min(AllScale_gp2), max(AllScale_gp2), 100);
    yy2 = polyval(p2, xx2);

    plot(xx2, yy2, ...
        'Color',CGROUP2*0.8, ...
        'LineWidth',2)

    % correlation
    [R2,P2] = corr(AllScale_gp2, AllField_gp2, ...
        'rows','complete');


    %% appearance

    xlabel('Grid scale')
    ylabel('Grid field width')

    legend({ ...
        sprintf('Ctrl  R=%.2f',R1), ...
        '', ...
        sprintf('GROUP2  R=%.2f',R2), ...
        ''}, ...
        'Box','off')

    box off

    title('All sessions pooled')


    if TrialSelection == 0
        s = strcat("EB2_Scatter_GridScaleVsWidth_All_260511.emf");
    else
        % s = strcat("EB2_Scatter_GridScaleVsWidth_ModuleTrial.emf");
    end

    xlim([20 80])
    xticks([20 40 60 80])
    % ylim([30 120])
    ylim([10 50])
    % yticks([30 40 60 80 100 120])
    yticks([10 30 50])
    exportgraphics(fig, s);



    %%
    %% =========================================================
    % Grid width vs scale (pooled + statistics)
    %% =========================================================

    pos2 = [300 500 160 160];
    fig = figure('Position',pos2);
    hold on

    % =========================
    % GROUP1
    % =========================
    Xn = AllField_gp1;
    Yn = AllScale_gp1;

    scatter(Xn, Yn, 10, ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',0.5)

    % linear fit
    p1 = polyfit(Xn, Yn, 1);
    xx1 = linspace(min(Xn), max(Xn), 100);
    yy1 = polyval(p1, xx1);
    plot(xx1, yy1, 'k', 'LineWidth',2)

    % correlation
    [R1,P1] = corr(Xn, Yn, 'rows','complete');


    % =========================
    % GROUP2
    % =========================
    Xc = AllField_gp2;
    Yc = AllScale_gp2;

    scatter(Xc, Yc, 10, ...
        'MarkerFaceColor','none', ...
        'MarkerEdgeColor',CGROUP2, ...
        'LineWidth',0.5)

    p2 = polyfit(Xc, Yc, 1);
    xx2 = linspace(min(Xc), max(Xc), 100);
    yy2 = polyval(p2, xx2);
    plot(xx2, yy2, 'Color',CGROUP2*0.8, 'LineWidth',2)

    [R2,P2] = corr(Xc, Yc, 'rows','complete');


    % =========================
    % Interaction model (KEY ANALYSIS)
    % =========================
    X = [Xn; Xc];
    Y = [Yn; Yc];
    Group = [zeros(size(Xn)); ones(size(Xc))]; % 0=GROUP1, 1=GROUP2

    Tbl = table(X, Y, categorical(Group), ...
        'VariableNames', {'X','Y','Group'});

    mdl = fitlm(Tbl, 'Y ~ X*Group');

    coef = mdl.Coefficients;

    % slope extraction
    beta_ctrl = coef.Estimate(2);        % X term
    beta_diff = coef.Estimate(4);        % interaction
    beta_gp2 = beta_ctrl + beta_diff;

    p_interaction = coef.pValue(4);


    % =========================
    % Print statistics
    % =========================
    fprintf('\n====================================\n')
    fprintf('Grid width vs scale analysis\n')
    fprintf('====================================\n')

    fprintf('GROUP1: R = %.4f, p = %.4f\n', R1, P1)
    fprintf('GROUP2:    R = %.4f, p = %.4f\n', R2, P2)

    fprintf('\n--- interaction model ---\n')
    fprintf('GROUP1 slope = %.6f\n', beta_ctrl)
    fprintf('GROUP2 slope    = %.6f\n', beta_gp2)
    fprintf('Interaction p = %.6f\n', p_interaction)

    fprintf('====================================\n')


    % =========================
    % Figure cosmetics
    % =========================
    ylabel('Grid scale')
    xlabel('Grid field width')

    legend({ ...
        sprintf('Ctrl  R=%.2f',R1), ...
        '', ...
        sprintf('GROUP2  R=%.2f',R2), ...
        ''}, ...
        'Box','off')

    box off

    xlim([10 50])
    xticks([10 30 50])
    ylim([20 80])
    yticks([20 40 60 80])

    xlabel('')
    ylabel('')

    fprintf('\n===== Linear regression =====\n')

    % linear regression model
    mdl_ctrl = fitlm(Xn, Yn);
    mdl_gp2 = fitlm(Xc, Yc);

    coef_ctrl = mdl_ctrl.Coefficients;
    coef_gp2 = mdl_gp2.Coefficients;

    % slope statistics
    slope_ctrl = coef_ctrl.Estimate(2);
    p_ctrl = coef_ctrl.pValue(2);
    t_ctrl = coef_ctrl.tStat(2);

    Slope_gp2 = coef_gp2.Estimate(2);
    p_gp2 = coef_gp2.pValue(2);
    t_gp2 = coef_gp2.tStat(2);

    % R squared
    R2_ctrl = mdl_ctrl.Rsquared.Ordinary;
    R2_gp2 = mdl_gp2.Rsquared.Ordinary;

    fprintf('GROUP1:\n')
    fprintf('slope = %.6f\n', slope_ctrl)
    fprintf('t = %.3f\n', t_ctrl)
    fprintf('p = %.6g\n', p_ctrl)
    fprintf('R^2 = %.4f\n', R2_ctrl)
    fprintf('R = %.4f, Pearson R = %.4f\n\n', sqrt(R2_ctrl), R1)

    fprintf('GROUP2:\n')
    fprintf('slope = %.6f\n', Slope_gp2)
    fprintf('t = %.3f\n', t_gp2)
    fprintf('p = %.6g\n', p_gp2)
    fprintf('R^2 = %.4f\n', R2_gp2)
    fprintf('R = %.4f, Pearson R = %.4f\n\n', sqrt(R2_gp2), R2)


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