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



%% GSratio

GridModNaoki =  load("..\08_GridModule\GridMod.mat");
A = GridModNaoki ;

GP2ratio = A.gp2_GSratio(:);
GP2ratio = GP2ratio(~isnan(GP2ratio));
GP1ratio = A.gp1_GSratio(:);
GP1ratio = GP1ratio(~isnan(GP1ratio));

[~, p_GP2] = ttest2(GP2ratio, GP1ratio);




%% ========== Compute Results of grid scale ==========

scale_mode = 'ratio'; % 'raw'
% scale_mode = 'raw'; % not scale ratio but raw scale [cm]
ResultsScale_Grp1  = compute_results_GridScale(GROUP1_T, A.gp1_GSpeak, A.gp1_GSratio, min_modulenum, mode, scale_mode);
ResultsScale_Grp2 = compute_results_GridScale(GROUP2_T, A.gp2_GSpeak, A.gp2_GSratio, min_modulenum, mode, scale_mode);




%% =================== Raw DV position vs Raw Grid Scale ===================

close all

% ===== Collect raw data =====
ScaleDV_Grp1_raw = [];
for Session_I = 1:height(ResultsScale_Grp1)

    % [raw grid scale, raw DV position, module ID]
    ScaleDV = [ ...
        ResultsScale_Grp1(Session_I).scale_ratio, ...
        ResultsScale_Grp1(Session_I).DV, ...
        ResultsScale_Grp1(Session_I).MI];

    ScaleDV_Grp1_raw = [ScaleDV_Grp1_raw; ScaleDV];

end


ScaleDV_Grp2_raw = [];
for Session_I = 1:height(ResultsScale_Grp2)

    ScaleDV = [ ...
        ResultsScale_Grp2(Session_I).scale_ratio, ...
        ResultsScale_Grp2(Session_I).DV, ...
        ResultsScale_Grp2(Session_I).MI];

    ScaleDV_Grp2_raw = [ScaleDV_Grp2_raw; ScaleDV];

end


% ===== Remove NaN =====
idx = ~any(isnan(ScaleDV_Grp1_raw(:,1:2)),2);
ScaleDV_Grp1_raw = ScaleDV_Grp1_raw(idx,:);

idx = ~any(isnan(ScaleDV_Grp2_raw(:,1:2)),2);
ScaleDV_Grp2_raw = ScaleDV_Grp2_raw(idx,:);


% Naoki's correction of DV position, based on plotColor_GridFieldVsLocation_Mod1norm.m
ScaleDV_Grp1_raw(:,2) = -2050 - ScaleDV_Grp1_raw(:,2);
ScaleDV_Grp2_raw(:,2) = -2050 - ScaleDV_Grp2_raw(:,2);

%% ===== Correlation =====
[r_Grp1, p_Grp1] = corr( ...
    ScaleDV_Grp1_raw(:,2), ...
    ScaleDV_Grp1_raw(:,1), ...
    'type','Pearson');

[r_Grp2, p_Grp2] = corr( ...
    ScaleDV_Grp2_raw(:,2), ...
    ScaleDV_Grp2_raw(:,1), ...
    'type','Pearson');


%% ===== Linear fit =====
mdl_Grp1 = fitlm( ...
    ScaleDV_Grp1_raw(:,2), ...
    ScaleDV_Grp1_raw(:,1));

mdl_Grp2 = fitlm( ...
    ScaleDV_Grp2_raw(:,2), ...
    ScaleDV_Grp2_raw(:,1));

slope_Grp1 = mdl_Grp1.Coefficients.Estimate(2);
slope_Grp2 = mdl_Grp2.Coefficients.Estimate(2);


%% =================== Scatter plot ===================

Cgroup1 = [0 0 0];
% Cgroup2 = [255 168 60]/255;
% Cgroup1 = [1 1 1];
Cgroup2 = [240 134 134]/255;

pos = [100 700 170 140];
fig = figure;
set(fig,'Position',pos)

hold on

% ===== Group1 =====
scatter( ...
    ScaleDV_Grp1_raw(:,2), ...
    ScaleDV_Grp1_raw(:,1), ...
    8, ...
    'MarkerEdgeColor',Cgroup1, ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.5)

% ===== Group2 =====
scatter( ...
    ScaleDV_Grp2_raw(:,2), ...
    ScaleDV_Grp2_raw(:,1), ...
    8, ...
    'MarkerEdgeColor',Cgroup2, ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.5)


% ===== Regression line =====
xfit = linspace( ...
    min([ScaleDV_Grp1_raw(:,2); ScaleDV_Grp2_raw(:,2)]), ...
    max([ScaleDV_Grp1_raw(:,2); ScaleDV_Grp2_raw(:,2)]), ...
    100);

yfit_Grp1 = predict(mdl_Grp1, xfit');
yfit_Grp2 = predict(mdl_Grp2, xfit');

plot(xfit, yfit_Grp1, ...
    'k', 'LineWidth',1)

plot(xfit, yfit_Grp2, ...
    'Color',Cgroup2, 'LineWidth',1)


box off

xlim([0 1000])


%% ===== Display statistics =====
fprintf('\n===== Group1 =====\n')
fprintf('r = %.3f\n', r_Grp1)
fprintf('p = %.5f\n', p_Grp1)
fprintf('slope = %.5f\n', slope_Grp1)

fprintf('\n===== Group2 =====\n')
fprintf('r = %.3f\n', r_Grp2)
fprintf('p = %.5f\n', p_Grp2)
fprintf('slope = %.5f\n', slope_Grp2)



%%

% ===== pooled data =====
X = [ ...
    ScaleDV_Grp1_raw(:,2);
    ScaleDV_Grp2_raw(:,2)];

Y = [ ...
    ScaleDV_Grp1_raw(:,1);
    ScaleDV_Grp2_raw(:,1)];

Group = [ ...
    zeros(size(ScaleDV_Grp1_raw,1),1);   % Group1
    ones(size(ScaleDV_Grp2_raw,1),1)];  % Group2

% 



%%
%% ===== Fixed intercept regression =====

% data from 5/9/2026
    AA(:,1) = readmatrix('..\08_GridModule\mingridscale.xlsx' ...
        , 'Sheet', 'gp1');

    AA(:,2) = readmatrix('..\08_GridModule\mingridscale.xlsx' ...
        , 'Sheet', 'gp2');

intercept_Grp1_all = AA(:,1);
intercept_Grp2_all = AA(:,2);
[~, p] = ttest2(intercept_Grp1_all, intercept_Grp2_all);
Intercept_Grp1  = mean(intercept_Grp1_all, 'omitmissing');
Intercept_Grp2  = mean(intercept_Grp2_all, 'omitmissing');

[p2,h,stats] = ranksum(intercept_Grp1_all, intercept_Grp2_all);

% ===== slope only =====
slope_Grp1 = regress( ...
    ScaleDV_Grp1_raw(:,1) - Intercept_Grp1, ...
    ScaleDV_Grp1_raw(:,2));

slope_Grp2 = regress( ...
    ScaleDV_Grp2_raw(:,1) - Intercept_Grp2, ...
    ScaleDV_Grp2_raw(:,2));


% ===== Regression line =====

xfit = linspace( ...
    min([ScaleDV_Grp1_raw(:,2); ScaleDV_Grp2_raw(:,2)]), ...
    max([ScaleDV_Grp1_raw(:,2); ScaleDV_Grp2_raw(:,2)]), ...
    100);

yfit_Grp1  = Intercept_Grp1  + slope_Grp1  * xfit;
yfit_Grp2 = Intercept_Grp2 + slope_Grp2 * xfit;

plot(xfit, yfit_Grp1, ...
    'k', 'LineWidth',1)

plot(xfit, yfit_Grp2, ...
    'Color',Cgroup2, 'LineWidth',1)


fprintf('\n===== Group1 =====\n')
fprintf('fixed intercept = %.1f\n', Intercept_Grp1)
fprintf('slope = %.5f\n', slope_Grp1)

fprintf('\n===== Group2 =====\n')
fprintf('fixed intercept = %.1f\n', Intercept_Grp2)
fprintf('slope = %.5f\n', slope_Grp2)




Yadj = Y;

Yadj(Group==0) = Yadj(Group==0) - Intercept_Grp1;
Yadj(Group==1) = Yadj(Group==1) - Intercept_Grp2;

Tbl = table(X, Yadj, categorical(Group));
Tbl.Properties.VariableNames = {'X','Yadj','Group'};

mdl_interaction = fitlm(Tbl, ...
    'Yadj ~ -1 + X + X:Group');

coefTable = mdl_interaction.Coefficients;
p_interaction = coefTable.pValue(2);

disp(p_interaction)





% =========================================================
% Plot fixed-intercept interaction model
% =========================================================

% ===== coefficients =====
coefTable = mdl_interaction.Coefficients;

% Group1 slope
beta_Grp1 = coefTable.Estimate(1);

% slope difference
beta_diff = coefTable.Estimate(2);

% Group2 slope
beta_Grp2 = beta_Grp1 + beta_diff;


% ===== x range =====
xfit = linspace( ...
    min(X), ...
    max(X), ...
    200);


% ===== fitted lines =====
yfit_Grp1 = Intercept_Grp1 + beta_Grp1 * xfit;

yfit_Grp2 = Intercept_Grp2 + beta_Grp2 * xfit;




% =========================================================
% Scatter + fitted line
% =========================================================

figure('Position',[100 700 170 140])
hold on

% ===== Group1 scatter =====
scatter( ...
    ScaleDV_Grp1_raw(:,2), ...
    ScaleDV_Grp1_raw(:,1), ...
    8, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.5)

% ===== Group2 scatter =====
scatter( ...
    ScaleDV_Grp2_raw(:,2), ...
    ScaleDV_Grp2_raw(:,1), ...
    8, ...
    'MarkerEdgeColor',Cgroup2, ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.5)


% ===== fitted lines =====
plot( ...
    xfit, ...
    yfit_Grp1, ...
    'k', ...
    'LineWidth',1.5)

plot( ...
    xfit, ...
    yfit_Grp2, ...
    'Color',Cgroup2, ...
    'LineWidth',1.5)


box off

xlabel('DV position')
ylabel('Grid scale ratio')

xlim([0 1000])


% ===== statistics =====

fprintf('\n====================================\n')
fprintf('Fixed-intercept slope comparison\n')
fprintf('Group1 slope = %.6f\n', beta_Grp1)
fprintf('Group2 slope    = %.6f\n', beta_Grp2)
fprintf('interaction p = %.6f\n', p_interaction)
fprintf('====================================\n')





%% Figure

% Cgroup2 = [240 134 134]/255;

fig = figure('Position',[100 700 100 140]);
hold on

% ===== Group1 scatter =====
scatter( ...
    ScaleDV_Grp1_raw(:,1), ...
    ScaleDV_Grp1_raw(:,2), ...
    6, ...
    'MarkerEdgeColor','k', ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.6)

% ===== Group2 scatter =====
scatter( ...
    ScaleDV_Grp2_raw(:,1), ...
    ScaleDV_Grp2_raw(:,2), ...
    6, ...
    'MarkerEdgeColor',Cgroup2*1.05, ...
    'MarkerFaceColor','none', ...
    'LineWidth',0.6)


% ===== fitted lines =====
plot(yfit_Grp1, xfit,  ...
    'k', 'LineWidth',2)

plot(yfit_Grp2, xfit, ...
    'Color',Cgroup2*0.9, ...
    'LineWidth',2)
% plot(yfit_Grp2, xfit, ...
%     'r', ...
%     'LineWidth',1.5)

box off

% ylabel('DV position')
% xlabel('Grid scale ratio')

ylim([0 1000])
xlim([0.5 2.5])
% xlim([0.5 2.0])
xticks([0.5 1.0 1.5 2.0 2.5])
yticks([0 200 400 600 800 1000])

ylabel('')
xlabel('')
xticklabels([])
yticklabels([])

set(gca, 'YDir','reverse')

H=gca;
H.LineWidth = 1.0;

s = strcat("DVposVsScale.emf");
% s = strcat("DVposVsScale_2.emf");
exportgraphics(fig, s);


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


        Gscale = OcT{m,d}.Grid_scale;
        Gmod   = OcT{m,d}.Grid_module;
        idx = (Gscale>0) & (Gmod>0);

        % ===== GSpeak =====
        GSpeak = GSpeak_all{m,d};
        if isempty(GSpeak) || isnan(GSratio_all(m,d))
            continue
        end
        module_peak_scale = GSpeak(:,1);
        if length(module_peak_scale) < min_modulenum
            continue
        end


        % ===== scale =====
        switch scale_mode
            case 'raw'
                [~, DV, MI] = func_grid_scale(OcT, m, d, module_peak_scale); 
                scale = Gscale(idx);
            case 'ratio'
                [scale_ratio, DV, MI] = func_grid_scale(OcT, m, d, module_peak_scale); 
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



