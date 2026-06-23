clear
close all

%% set parameters

Base = 'H:\experiments H drive\BayesResults_260317_BorderCenter';

celltype_str = {"place cells", "all reliable cells", "all cells"};

OutDir = 'H:\experiments H drive\Figs_260504_BorderCenter_AllDecoder';
[status,msg,msgID] = mkdir(fullfile(OutDir));


%% load data

CellN = 5:1:60;

N = length(CellN);

Results = cell(N,1);
for i = 1:N

    try
        cname = sprintf('%dcells', CellN(i));
        fprintf('load %s data\n', cname)
        Results{i} = load(fullfile(Base, cname, 'bayesian_decoding_results.mat'));
    catch
    end
end

SessionList = Results{1}.session_list;
Nsession =  size(SessionList, 2);



%% ====== decoding error (All decoder + distance-based BC) ======
close all

fig_cell_all      = cell(3,3);
fig_diff_cell_all = cell(3,3);
fig_Z_cell_all    = cell(3,3);

Match_M_all       = cell(3,3);

arena_x = 100;
arena_y = 100;
th = 20; % cm

BC_str = {'All','Center','Border'};

for BC_ind = 1:3 %1: All, 2: Center, 3: Border
    
    Match_M_celltype = cell(1,3);
    Match_M_shuf_celltype = cell(1,3);
    Match_Z_shuf_celltype = cell(1,3);

    for celltype_i_temp = 1 %1:3
        
        Match_M        = nan(Nsession, length(Results));
        Match_M_shuf   = nan(Nsession, length(Results));
        Match_Z_shuf   = nan(Nsession, length(Results));

        for i = 1:length(Results)

            if isempty(Results{i})
                continue
            end

            tmp = Results{i}.all_decoding_results{celltype_i_temp,1};

            vals_real = nan(length(tmp),1);
            vals_shuf = nan(length(tmp),1);
            vals_std  = nan(length(tmp),1);

            for j = 1:length(tmp)

                if ~isstruct(tmp{j})
                    continue
                end

                try
                    % ===== position =====
                    pos = tmp{j}.observed_xy;
                    x = pos(:,1); y = pos(:,2);
                    dist = min([x, arena_x-x, y, arena_y-y], [], 2);


                    % ===== mask =====
                    if BC_ind == 1
                        mask = true(size(dist));  % All
                    elseif BC_ind == 2
                        mask = dist > th;         % Center
                    elseif BC_ind == 3
                        mask = dist <= th;        % Border
                    end


                    A = tmp{j}.matched_per_frame_errors;
                    nSub = size(A,2);
                    sub_vals = nan(nSub,1);
                    for s = 1:nSub
                        err_s = A(:,s); % each subsample
                        err_s = err_s(mask);
                        sub_vals(s) = median(err_s, 'omitnan');
                    end
                    vals_real(j) = mean(sub_vals, 'omitnan');
       

                    % ===== SHUFFLE=====
                    if isfield(tmp{j}, 'shuf_matched_per_frame_errors_shuf')

                        A_shuf = tmp{j}.shuf_matched_per_frame_errors_shuf;

                        if iscell(A_shuf)
                            A_shuf = cell2mat(cellfun(@(x) x(:), A_shuf, 'UniformOutput', false));
                        end

             
                        nSub = size(A_shuf,2);
                        shuf_sub = nan(nSub,1);
                        for s = 1:nSub
                            err_s = A_shuf(:,s); % each subsample
                            err_s = err_s(mask);
                            shuf_sub(s) = median(err_s, 'omitnan');
                        end

                        vals_shuf(j) = mean(shuf_sub, 'omitnan');
                        vals_std(j)  = std(shuf_sub, 'omitnan');
                    end

                catch
                end
            end

            Match_M(:, i)      = vals_real;
            Match_M_shuf(:, i) = vals_shuf;

            % ===== Z-score =====
            Match_Z_shuf(:, i) = (vals_real - vals_shuf) ./ vals_std;

        end

        Match_M_celltype{celltype_i_temp}        = Match_M;
        Match_M_shuf_celltype{celltype_i_temp}   = Match_M_shuf;
        Match_Z_shuf_celltype{celltype_i_temp}   = Match_Z_shuf;

        Match_M_all{BC_ind, celltype_i_temp} = Match_M;
    end

    % ===== Plot =====
    for celltype_i_temp = 1 %1:3

        Cont_ind = find(strcmp({SessionList.group_label}, 'WildType'));
        Casp_ind = find(strcmp({SessionList.group_label}, 'IslandKilled'));

        Match_M = Match_M_celltype{celltype_i_temp};
        Control = Match_M(Cont_ind,:);
        Casp    = Match_M(Casp_ind,:);

        Match_M_shuf = Match_M_shuf_celltype{celltype_i_temp};
        shuf_Control = Match_M_shuf(Cont_ind,:);
        shuf_Casp    = Match_M_shuf(Casp_ind,:);

        Match_Z = Match_Z_shuf_celltype{celltype_i_temp};
        Control_Z = Match_Z(Cont_ind,:);
        Casp_Z    = Match_Z(Casp_ind,:);


        fprintf(strcat("================ DecodingError, ", BC_str{BC_ind}, "\n")  );
        [fig, lme, fig_diff, fig_Z] = func_PlotFitting( ...
            Control, Casp, CellN, celltype_i_temp, celltype_str, ...
            shuf_Control, shuf_Casp, SessionList, Control_Z, Casp_Z, BC_ind);

        fig_cell_all{BC_ind, celltype_i_temp}      = fig;
        fig_diff_cell_all{BC_ind, celltype_i_temp} = fig_diff;
        fig_Z_cell_all{BC_ind, celltype_i_temp}    = fig_Z;

    end
end


% ===== save figs =====
for BC_ind = 1:3
    for i = 1 %1:3

        s = strcat("DecodingError_", BC_str{BC_ind}, "_", celltype_str{i}, '.pdf');
        exportgraphics(fig_cell_all{BC_ind,i}, fullfile(OutDir, s));

        s = strcat("DecodingErrorZ_", BC_str{BC_ind}, "_", celltype_str{i}, '.pdf');
        exportgraphics(fig_Z_cell_all{BC_ind,i}, fullfile(OutDir, s));

    end
end

close all




%% =========== Accuracy (All decoder + distance-based BC) ==============
close all

Acc_Bin_list = [15 30];
BC_str = {'All','Center','Border'};

Acc_Ratio_all = cell(3,3,length(Acc_Bin_list)); % BC × celltype × AccBin
fig_acc_all   = cell(3,3,length(Acc_Bin_list));

% ===== arena & threshold =====
arena_x = 100;
arena_y = 100;
th = 20; % cm (boundary definition)

for b = 1 %1:length(Acc_Bin_list)
    Acc_Bin = Acc_Bin_list(b);

    for BC_ind = 1:3 % 1: All, 2: Center, 3: Border
        Acc_Ratio_celltype = cell(1,3);

        for celltype_i_temp = 1:3
            Acc_Ratio = nan(Nsession, length(Results));

            for i = 1:length(Results)

                if isempty(Results{i})
                    continue
                end

                tmp = Results{i}.all_decoding_results{celltype_i_temp,1};
                vals = nan(length(tmp),1);

                for j = 1:length(tmp)

                    if ~isstruct(tmp{j})
                        continue
                    end

                    try
                        %%
                        % ===== position =====
                        pos = tmp{j}.observed_xy;
                        x = pos(:,1);
                        y = pos(:,2);

                        dist = min([x, arena_x-x, y, arena_y-y], [], 2);

                        % ===== mask =====
                        if BC_ind == 1
                            mask = true(size(dist));  % All
                        elseif BC_ind == 2
                            mask = dist > th;         % Center
                        elseif BC_ind == 3
                            mask = dist <= th;        % Border
                        end

                        err = tmp{j}.matched_per_frame_errors;

                        if BC_ind == 1
                            err_masked = err(~isnan(err));
                        else
                            err_masked = err(mask,:);
                            err_masked = err_masked(~isnan(err_masked));
                        end

                        % Accuracy
                        AA = sum(err_masked(:) < Acc_Bin) / length(err_masked(:));
                        vals(j) = AA;

                    catch
                    end
                end

                Acc_Ratio(1:length(vals), i) = vals;
            end

            Acc_Ratio_celltype{celltype_i_temp} = Acc_Ratio;
            Acc_Ratio_all{BC_ind, celltype_i_temp, b} = Acc_Ratio;
        end


        % ===== Plot =====
        fprintf( strcat('\n============', BC_str{BC_ind}, '\n'))
        for celltype_i_temp = 1 %1:3

            fig = func_Accuracy(Acc_Ratio_celltype, celltype_i_temp, SessionList, CellN, celltype_str, b);
            fig_acc_all{BC_ind, celltype_i_temp, b} = fig;
        end

    end
end

%
% save
for b = 1 %1:length(Acc_Bin_list)
    Acc_Bin = Acc_Bin_list(b);
    for BC_ind = 1:3
        for i = 1 %1:3

            s = sprintf("Accuracy_%s_%dcm_%s.pdf", ...
                BC_str{BC_ind}, Acc_Bin, celltype_str{i});

            exportgraphics(fig_acc_all{BC_ind, i, b}, fullfile(OutDir, s));

        end
    end
end




%% functions
%% plot fitting results
function [fig, lme, fig_diff, fig_Z] = func_PlotFitting(Control, Casp, CellN, celltype_i_temp, celltype_str, shuf_Control, shuf_Casp, SessionList, Control_Z, Casp_Z, BC_ind)

Ylim = [0 65];

%% plot results
celltype = celltype_i_temp;


Cont_color = [1 1 1]*0.2;
Casp_color =[255 168 60]/255;

M_Cont = mean(Control, 'omitmissing');
M_Casp = mean(Casp, 'omitmissing');

SEM_Cont = std(Control, 0, 'omitmissing') ./ sqrt(sum(~isnan(Control)));
SEM_Casp = std(Casp, 0, 'omitmissing') ./ sqrt(sum(~isnan(Casp)));


% get(gcf,'Position')
pos = [65   678   810   275];
fig = figure('Position', pos);
% % fig = figure;
% subplot(1,2,1)
hold on

% --- scatter (each trial) ---
for i = 1:length(CellN)
    jitter = 0.5;
    scatter(CellN(i) - jitter*randn(size(Control(:,i))), Control(:,i), ...
        10, Cont_color, 'filled', 'MarkerFaceAlpha', 0.3)
    scatter(CellN(i) + jitter*randn(size(Casp(:,i))), Casp(:,i), ...
        10, Casp_color, 'filled', 'MarkerFaceAlpha', 0.3)
end


plot(CellN, M_Cont, 'Color', Cont_color, 'LineWidth', 3)
plot(CellN, M_Casp, 'Color', Casp_color, 'LineWidth', 3)

errorbar(CellN, M_Cont, SEM_Cont, ...
    'Color', Cont_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)
errorbar(CellN, M_Casp, SEM_Casp, ...
    'Color', Casp_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)

ylim(Ylim)
xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('Mean decoding error [cm]')
box off
set(gca,'FontSize',10)


% plot shuffle data
if isempty(shuf_Control) == 0 && isempty(shuf_Casp) == 0
    shuf_Cont_color = [1 1 1] * 0.8;
    shuf_Casp_color = [1, 0.85, 0.7];

    shuf_M_Cont = mean(shuf_Control, 'omitmissing');
    shuf_M_Casp = mean(shuf_Casp, 'omitmissing');

    shuf_SEM_Cont = std(shuf_Control, 0, 'omitmissing') ./ sqrt(sum(~isnan(shuf_Control)));
    shuf_SEM_Casp = std(shuf_Casp, 0, 'omitmissing') ./ sqrt(sum(~isnan(shuf_Casp)));

    for i = 1:length(CellN)
        jitter = 0.5;
        scatter(CellN(i) - jitter*randn(size(shuf_Control(:,i))), ...
            shuf_Control(:,i), ...
            10, shuf_Cont_color, 'filled', ...
            'MarkerFaceAlpha', 0.3)
        scatter(CellN(i) + jitter*randn(size(shuf_Casp(:,i))), ...
            shuf_Casp(:,i), ...
            10, shuf_Casp_color, 'filled', ...
            'MarkerFaceAlpha', 0.3)
    end

    e1 = errorbar(CellN, shuf_M_Cont, shuf_SEM_Cont, ...
        'Color', shuf_Cont_color + 0, ...
        'LineStyle', 'none', ...
        'LineWidth', 1.5);
    e2 = errorbar(CellN, shuf_M_Casp, shuf_SEM_Casp, ...
        'Color', shuf_Casp_color + 0.0, ...
        'LineStyle', 'none', ...
        'LineWidth', 1.5);

    p1 = plot(CellN, shuf_M_Cont, 'Color', shuf_Cont_color, 'LineWidth', 3);
    p2 = plot(CellN, shuf_M_Casp, 'Color', shuf_Casp_color, 'LineWidth', 3);
    p1.Color(4) = 0.4;
    p2.Color(4) = 0.5;


    [p_Control,h,stats] = ranksum(Control(:),shuf_Control(:));
    [p_Casp,h,stats] = ranksum(Casp(:),shuf_Casp(:));
    title(sprintf('Control vs. Control_shuffle: p = %.4f\nCasp vs. Casp_shuffle = %.4f', p_Control, p_Casp), 'Interpreter','none');

end

%% create fitting model

% group matrices
% Control (animals x N)
% Casp (animals x N)

[nA_ctrl, nN] = size(Control);
[nA_casp, nN2] = size(Casp);

InvSqrtN = 1 ./ sqrt(CellN);

% table
Error = [];
InvSqrtN_all = [];

Group = strings(0);
Session = strings(0);

Animal = strings(0);
RecDay = strings(0);

% Control
for a = 1:nA_ctrl
    for n = 1:nN
        val = Control(a,n);
        if ~isnan(val)
            Error(end+1,1) = val;
            InvSqrtN_all(end+1,1) = InvSqrtN(n);
            Group(end+1,1) = "Control";
            Session(end+1,1) = "C" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a);
        end
    end
end

% Casp
for a = 1:nA_casp
    for n = 1:nN2
        val = Casp(a,n);
        if ~isnan(val)
            Error(end+1,1) = val;
            InvSqrtN_all(end+1,1) = InvSqrtN(n);
            Group(end+1,1) = "Casp";
            Session(end+1,1) = "K" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a + nA_ctrl);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a + nA_ctrl);
        end
    end
end

% remove invalid rows
valid_idx = ...
    ~isnan(Error) & ...
    ~isnan(InvSqrtN_all) & ...
    ~isinf(Error) & ...
    ~isinf(InvSqrtN_all);

Error = Error(valid_idx);
InvSqrtN_all = InvSqrtN_all(valid_idx);
Group = Group(valid_idx);
% Session = Session(valid_idx);
Animal = Animal(valid_idx);

T = table(Error,InvSqrtN_all,categorical(Group),categorical(Animal), categorical(RecDay), categorical(Session), ...
    'VariableNames',{'Error','InvSqrtN','Group','Animal', 'RecDay', 'Session'});
T = rmmissing(T);

lme = fitlme(T,'Error ~ InvSqrtN*Group + (InvSqrtN|Session)');
fprintf('BC_ind == %d, CellType: %s\n', BC_ind, celltype_str{celltype})
anova(lme,'DFMethod','Satterthwaite')


%% print statistics
aov = anova(lme,'DFMethod','Satterthwaite');

fprintf('\n===== Decoding error (cm) =====\n')
fprintf('\n===== Mixed-effects model =====\n')

% split table
T_ctrl = T(T.Group == "Control", :);
T_casp = T(T.Group == "Casp", :);

% Control
n_obs_ctrl = height(T_ctrl);
n_session_ctrl = numel(unique(T_ctrl.Session));
n_animal_ctrl = numel(unique(T_ctrl.Animal));

% Casp
n_obs_casp = height(T_casp);
n_session_casp = numel(unique(T_casp.Session));
n_animal_casp = numel(unique(T_casp.Animal));

% print
fprintf('\n===== Sample size =====\n')
fprintf('Control: n = %d observations, %d sessions, %d animals\n', ...
    n_obs_ctrl, n_session_ctrl, n_animal_ctrl);
fprintf('Casp:    n = %d observations, %d sessions, %d animals\n', ...
    n_obs_casp, n_session_casp, n_animal_casp);


% stats
for i = 1:height(aov)
    term = string(aov.Term(i));
    F = aov.FStat(i);
    df1 = aov.DF1(i);
    df2 = aov.DF2(i);
    p = aov.pValue(i);
    if p < 0.0001
        pstr = "p < 0.0001";
    else
        pstr = sprintf("p = %.4g", p);
    end
    fprintf('%s: %s, F (%.3f, %.3f) = %.3f\n', ...
        term, pstr, df1, df2, F);
end

%% visualization

% set(0, 'currentfigure', fig);

Nfit = linspace(min(CellN), max(CellN), 200)';
InvFit = 1 ./ sqrt(Nfit);

Tpred_Control = table(InvFit, ...
    categorical(repmat("Control",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'InvSqrtN','Group','Session'}); %All Session is 1 -> predict only population-level curve

Tpred_Casp = table(InvFit, ...
    categorical(repmat("Casp",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'InvSqrtN','Group','Session'});

yfit_Control = predict(lme,Tpred_Control);
yfit_Casp = predict(lme,Tpred_Casp);


% figure
subplot(1,2,2)
hold on


% individual points
MarkerSize = 4;
Alpha = 0.2;
for i = 1:length(CellN)
    % jitter = 0.6;
    scatter(CellN(i) - jitter*randn(size(Control(:,i))), ...
        Control(:,i), MarkerSize, Cont_color, ...
        'filled','MarkerFaceAlpha', Alpha)
    scatter(CellN(i) + jitter*randn(size(Casp(:,i))), ...
        Casp(:,i), MarkerSize, Casp_color, ...
        'filled','MarkerFaceAlpha', Alpha)
end


% model fit
plot(Nfit,yfit_Control,'Color',Cont_color,'LineWidth',2)
plot(Nfit,yfit_Casp,'Color',Casp_color,'LineWidth',2)


ylim(Ylim)
box off
set(gca,'FontSize',10)


anova_tbl = anova(lme,'DFMethod','Satterthwaite'); % table
p_Group = anova_tbl.pValue(strcmp(anova_tbl.Term,'Group'));
p_Int = anova_tbl.pValue(strcmp(anova_tbl.Term,'InvSqrtN:Group'));
% title(sprintf('Population curve\np_Group=%.3f, p_Interaction=%.3f', p_Group, p_Int), 'Interpreter','none');


xticks([0 30 60])
yticks([0 30 60])



%% plot shuffle data
if isempty(shuf_Control) == 0 && isempty(shuf_Casp) == 0
    % create fitting model

    [nA_ctrl, nN] = size(shuf_Control);
    [nA_casp, nN2] = size(shuf_Casp);

    InvSqrtN = 1 ./ sqrt(CellN);

    % table
    Error = [];
    InvSqrtN_all = [];
    Group = strings(0);
    Session = strings(0);

    % Control
    for a = 1:nA_ctrl
        for n = 1:nN
            val = shuf_Control(a,n);
            if ~isnan(val)
                Error(end+1,1) = val;
                InvSqrtN_all(end+1,1) = InvSqrtN(n);
                Group(end+1,1) = "shuf_Control";
                Session(end+1,1) = "C" + a;
            end
        end
    end

    % Casp
    for a = 1:nA_casp
        for n = 1:nN2
            val = shuf_Casp(a,n);
            if ~isnan(val)
                Error(end+1,1) = val;
                InvSqrtN_all(end+1,1) = InvSqrtN(n);
                Group(end+1,1) = "shuf_Casp";
                Session(end+1,1) = "K" + a;
            end
        end
    end

    shuf_CellAll = [Error, InvSqrtN_all, Group, Session];

    % remove invalid rows
    valid_idx = ...
        ~isnan(Error) & ...
        ~isnan(InvSqrtN_all) & ...
        ~isinf(Error) & ...
        ~isinf(InvSqrtN_all);

    Error = Error(valid_idx);
    InvSqrtN_all = InvSqrtN_all(valid_idx);
    Group = Group(valid_idx);
    Session = Session(valid_idx);

    T = table(Error,InvSqrtN_all,categorical(Group),categorical(Session),...
        'VariableNames',{'Error','InvSqrtN','Group','Session'});
    T = rmmissing(T);
    shuf_lme = fitlme(T,'Error ~ InvSqrtN*Group + (InvSqrtN|Session)');


    Nfit = linspace(min(CellN), max(CellN), 200)';
    InvFit = 1 ./ sqrt(Nfit);

    Tpred_Control = table(InvFit, ...
        categorical(repmat("shuf_Control",length(Nfit),1)), ...
        categorical(repmat("C1",length(Nfit),1)), ...
        'VariableNames',{'InvSqrtN','Group','Session'});

    Tpred_Casp = table(InvFit, ...
        categorical(repmat("shuf_Casp",length(Nfit),1)), ...
        categorical(repmat("K1",length(Nfit),1)), ...
        'VariableNames',{'InvSqrtN','Group','Session'});

    yfit_Control = predict(shuf_lme,Tpred_Control);
    yfit_Casp = predict(shuf_lme,Tpred_Casp);


    % figure
    subplot(1,2,2)
    hold on

    % individual points
    for i = 1:length(CellN)
        % jitter = 0.6;
        s1 = scatter(CellN(i) - jitter*randn(size(shuf_Control(:,i))), ...
            shuf_Control(:,i), MarkerSize, shuf_Cont_color, ...
            'filled','MarkerFaceAlpha',Alpha + 0.0);
        s2 = scatter(CellN(i) + jitter*randn(size(shuf_Casp(:,i))), ...
            shuf_Casp(:,i), MarkerSize, shuf_Casp_color, ...
            'filled','MarkerFaceAlpha',Alpha + 0.0);
    end


    % model fit
    h1 = plot(Nfit,yfit_Control,'Color',shuf_Cont_color,'LineWidth',2);
    h2 = plot(Nfit,yfit_Casp,'Color',shuf_Casp_color,'LineWidth',2);
    h1.Color(4) = 0.4;
    h2.Color(4) = 0.5;

    uistack(h1, 'bottom');
    uistack(h2, 'bottom');
    uistack(s1, 'bottom');
    uistack(s2, 'bottom');

    % xlabel('Number of cells (N)')
    xlabel(strcat("Number of ", celltype_str{celltype}') )
    ylabel('Mean decoding error (cm)')

    % ylim([0 Ymax])
    box off
    set(gca,'FontSize',10)
end




%% Difference with shuffle
%% plot results
Ylim = [-50 5];

celltype = celltype_i_temp;

Cont_diff = Control - shuf_Control;
Casp_diff = Casp - shuf_Casp;

M_Cont = mean(Cont_diff, 'omitmissing');
M_Casp = mean(Casp_diff, 'omitmissing');

SEM_Cont = std(Cont_diff, 0, 'omitmissing') ./ sqrt(sum(~isnan(Cont_diff)));
SEM_Casp = std(Casp_diff, 0, 'omitmissing') ./ sqrt(sum(~isnan(Casp_diff)));

fig_diff = figure('Position', pos);
% fig_diff = figure;
% subplot(1,2,1)
hold on

% --- scatter  ---
for i = 1:length(CellN)
    jitter = 0.5;
    scatter(CellN(i) - jitter*0.3*randn(size(Cont_diff(:,i))), Cont_diff(:,i), ...
        MarkerSize, Cont_color, 'filled', 'MarkerFaceAlpha', Alpha)
    scatter(CellN(i) + jitter*0.3*randn(size(Casp_diff(:,i))), Casp_diff(:,i), ...
        MarkerSize, Casp_color, 'filled', 'MarkerFaceAlpha', Alpha)
end

plot(CellN, M_Cont, 'Color', Cont_color, 'LineWidth', 3)
plot(CellN, M_Casp, 'Color', Casp_color, 'LineWidth', 3)

errorbar(CellN, M_Cont, SEM_Cont, ...
    'Color', Cont_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)

errorbar(CellN, M_Casp, SEM_Casp, ...
    'Color', Casp_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)

ylim(Ylim)
xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('Mean decoding error improvement [cm]')
box off
set(gca,'FontSize',10)


%% create fitting model

% group matrices
% Control (animals x N)
% Casp (animals x N)

[nA_ctrl, nN] = size(Cont_diff);
[nA_casp, nN2] = size(Casp_diff);

InvSqrtN = 1 ./ sqrt(CellN);

% table
Error = [];
InvSqrtN_all = [];
Group = strings(0);
Session = strings(0);

Animal = strings(0);
RecDay = strings(0);

% Control
for a = 1:nA_ctrl
    for n = 1:nN
        val = Cont_diff(a,n);
        if ~isnan(val)
            Error(end+1,1) = val;
            InvSqrtN_all(end+1,1) = InvSqrtN(n);
            Group(end+1,1) = "Control";
            Session(end+1,1) = "C" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a);
        end
    end
end

% Casp
for a = 1:nA_casp
    for n = 1:nN2
        val = Casp_diff(a,n);
        if ~isnan(val)
            Error(end+1,1) = val;
            InvSqrtN_all(end+1,1) = InvSqrtN(n);
            Group(end+1,1) = "Casp";
            Session(end+1,1) = "K" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a + nA_ctrl);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a + nA_ctrl);
        end
    end
end


% remove invalid rows
valid_idx = ...
    ~isnan(Error) & ...
    ~isnan(InvSqrtN_all) & ...
    ~isinf(Error) & ...
    ~isinf(InvSqrtN_all);

Error = Error(valid_idx);
InvSqrtN_all = InvSqrtN_all(valid_idx);
Group = Group(valid_idx);
% Session = Session(valid_idx);
Animal = Animal(valid_idx);


T = table(Error,InvSqrtN_all,categorical(Group),categorical(Animal), categorical(RecDay), categorical(Session), ...
    'VariableNames',{'Error','InvSqrtN','Group','Animal', 'RecDay', 'Session'});
T = rmmissing(T);


lme = fitlme(T,'Error ~ InvSqrtN*Group + (InvSqrtN|Session)');
anova(lme,'DFMethod','Satterthwaite')


%% visualization

Nfit = linspace(min(CellN), max(CellN), 200)';
InvFit = 1 ./ sqrt(Nfit);

Tpred_Control = table(InvFit, ...
    categorical(repmat("Control",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'InvSqrtN','Group','Session'}); %All Session is 1 -> predict only population-level curve

Tpred_Casp = table(InvFit, ...
    categorical(repmat("Casp",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'InvSqrtN','Group','Session'});

yfit_Control = predict(lme,Tpred_Control);
yfit_Casp = predict(lme,Tpred_Casp);


% figure
subplot(1,2,2)
hold on


% individual points
for i = 1:length(CellN)
    jitter = 0.6;
    scatter(CellN(i) - jitter*randn(size(Cont_diff(:,i))), ...
        Cont_diff(:,i), MarkerSize, Cont_color, ...
        'filled','MarkerFaceAlpha',Alpha)
    scatter(CellN(i) + jitter*randn(size(Casp_diff(:,i))), ...
        Casp_diff(:,i), MarkerSize, Casp_color, ...
        'filled','MarkerFaceAlpha',Alpha)
end


% model fit
plot(Nfit,yfit_Control,'Color',Cont_color,'LineWidth',3)
plot(Nfit,yfit_Casp,'Color',Casp_color,'LineWidth',3)

xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('Mean decoding error improvement [cm]')

% ylim([-30 5])
ylim(Ylim)
box off
set(gca,'FontSize',10)


anova_tbl = anova(lme,'DFMethod','Satterthwaite'); % 型は table
p_Group = anova_tbl.pValue(strcmp(anova_tbl.Term,'Group'));
p_Int = anova_tbl.pValue(strcmp(anova_tbl.Term,'InvSqrtN:Group'));
title(sprintf('Population curve\np_Group=%.3f, p_Interaction=%.3f', p_Group, p_Int), 'Interpreter','none');




%% ========== Z-score

Ylim = [-60 12];

% plot results
celltype = celltype_i_temp;

Cont_Z = Control_Z;

Z_Cont = mean(Cont_Z, 'omitmissing');
Z_Casp = mean(Casp_Z, 'omitmissing');

SEM_Cont = std(Cont_Z, 0, 'omitmissing') ./ sqrt(sum(~isnan(Cont_Z)));
SEM_Casp = std(Casp_Z, 0, 'omitmissing') ./ sqrt(sum(~isnan(Casp_Z)));


fig_Z = figure('Position', pos);
hold on

% --- scatter ---
for i = 1:length(CellN)
    jitter = 0.5;
    scatter(CellN(i) - jitter*0.3*randn(size(Cont_Z(:,i))), Cont_Z(:,i), ...
        MarkerSize, Cont_color, 'filled', 'MarkerFaceAlpha', Alpha)
    scatter(CellN(i) + jitter*0.3*randn(size(Casp_Z(:,i))), Casp_Z(:,i), ...
        MarkerSize, Casp_color, 'filled', 'MarkerFaceAlpha', Alpha)
end

plot(CellN, Z_Cont, 'Color', Cont_color, 'LineWidth', 2)
plot(CellN, Z_Casp, 'Color', Casp_color, 'LineWidth', 2)

errorbar(CellN, Z_Cont, SEM_Cont, ...
    'Color', Cont_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)
errorbar(CellN, Z_Casp, SEM_Casp, ...
    'Color', Casp_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)

ylim(Ylim)
xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('decoding error (z)')
box off
set(gca,'FontSize',10)


%% create fitting model

% group matrices
% Control (animals x N)
% Casp (animals x N)

[nA_ctrl, nN] = size(Cont_Z);
[nA_casp, nN2] = size(Casp_Z);

InvSqrtN = 1 ./ sqrt(CellN);

% table
Error = [];
InvSqrtN_all = [];

Group = strings(0);
Session = strings(0);

Animal = strings(0);
RecDay = strings(0);

% Control
for a = 1:nA_ctrl
    for n = 1:nN
        val = Cont_Z(a,n);
        if ~isnan(val)
            Error(end+1,1) = val;
            InvSqrtN_all(end+1,1) = InvSqrtN(n);
            Group(end+1,1) = "Control";
            Session(end+1,1) = "C" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a);
        end
    end
end

% Casp
for a = 1:nA_casp
    for n = 1:nN2
        val = Casp_Z(a,n);
        if ~isnan(val)
            Error(end+1,1) = val;
            InvSqrtN_all(end+1,1) = InvSqrtN(n);
            Group(end+1,1) = "Casp";
            Session(end+1,1) = "K" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a + nA_ctrl);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a + nA_ctrl);
        end
    end
end

% remove invalid rows
valid_idx = ...
    ~isnan(Error) & ...
    ~isnan(InvSqrtN_all) & ...
    ~isinf(Error) & ...
    ~isinf(InvSqrtN_all);

Error = Error(valid_idx);
InvSqrtN_all = InvSqrtN_all(valid_idx);
Group = Group(valid_idx);
% % Session = Session(valid_idx);
Animal = Animal(valid_idx);
RecDay = RecDay(valid_idx);
Session = Session(valid_idx);

T = table(Error,InvSqrtN_all,categorical(Group),categorical(Animal), categorical(RecDay), categorical(Session), ...
    'VariableNames',{'Error','InvSqrtN','Group','Animal', 'RecDay', 'Session'});

T = rmmissing(T);

lme = fitlme(T,'Error ~ InvSqrtN*Group + (InvSqrtN|Session)');
anova(lme,'DFMethod','Satterthwaite')


%%
%% print statistics
aov = anova(lme,'DFMethod','Satterthwaite');

fprintf('\n===== Decoding error (z-scored) =====\n')
fprintf('\n===== Mixed-effects model =====\n')


% split table
T_ctrl = T(T.Group == "Control", :);
T_casp = T(T.Group == "Casp", :);

% Control
n_obs_ctrl = height(T_ctrl);
n_session_ctrl = numel(unique(T_ctrl.Session));
n_animal_ctrl = numel(unique(T_ctrl.Animal));

% Casp
n_obs_casp = height(T_casp);
n_session_casp = numel(unique(T_casp.Session));
n_animal_casp = numel(unique(T_casp.Animal));

% print
fprintf('\n===== Sample size =====\n')
fprintf('Control: n = %d observations, %d sessions, %d animals\n', ...
    n_obs_ctrl, n_session_ctrl, n_animal_ctrl);
fprintf('Casp:    n = %d observations, %d sessions, %d animals\n', ...
    n_obs_casp, n_session_casp, n_animal_casp);


% stats
for i = 1:height(aov)
    term = string(aov.Term(i));
    F = aov.FStat(i);
    df1 = aov.DF1(i);
    df2 = aov.DF2(i);
    p = aov.pValue(i);

    if p < 0.0001
        pstr = "p < 0.0001";
    else
        pstr = sprintf("p = %.4g", p);
    end

    fprintf('%s: %s, F (%.3f, %.3f) = %.3f\n', ...
        term, pstr, df1, df2, F);
end

%% visualization

% set(0, 'currentfigure', fig);
Nfit = linspace(min(CellN), max(CellN), 200)';
InvFit = 1 ./ sqrt(Nfit);

Tpred_Control = table(InvFit, ...
    categorical(repmat("Control",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'InvSqrtN','Group','Session'}); %All Session is 1 -> predict only population-level curve

Tpred_Casp = table(InvFit, ...
    categorical(repmat("Casp",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'InvSqrtN','Group','Session'});

yfit_Control = predict(lme,Tpred_Control);
yfit_Casp = predict(lme,Tpred_Casp);


% figure
subplot(1,2,2)
hold on


% individual points
for i = 1:length(CellN)
    jitter = 0.6;
    scatter(CellN(i) - jitter*randn(size(Cont_Z(:,i))), ...
        Cont_Z(:,i), MarkerSize, Cont_color, ...
        'filled','MarkerFaceAlpha',Alpha)
    scatter(CellN(i) + jitter*randn(size(Casp_Z(:,i))), ...
        Casp_Z(:,i), MarkerSize, Casp_color, ...
        'filled','MarkerFaceAlpha',Alpha)
end


% model fit
plot(Nfit,yfit_Control,'Color',Cont_color,'LineWidth',2)
plot(Nfit,yfit_Casp,'Color',Casp_color,'LineWidth',2)

% ylim([-105 12])
ylim(Ylim)
xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('decoding error (z)')

box off
set(gca,'FontSize',10)


anova_tbl = anova(lme,'DFMethod','Satterthwaite'); % 型は table
p_Group = anova_tbl.pValue(strcmp(anova_tbl.Term,'Group'));
p_Int = anova_tbl.pValue(strcmp(anova_tbl.Term,'InvSqrtN:Group'));
% title(sprintf('Population curve\np_Group=%.3f, p_Interaction=%.3f', p_Group, p_Int), 'Interpreter','none');


if BC_ind == 1
    yticks([-100 -50 0])
elseif BC_ind == 2
    yticks([-60 -30 0])
elseif BC_ind == 3
    yticks([-100 -50 0])
end

xticks([0 30 60])
yticks([-60 -30 0])

end




%% plot Accuracy
function fig = func_Accuracy(Acc_Ratio_celltype, celltype_i_temp, SessionList, CellN, celltype_str, b)
%%
Cont_ind = find(strcmp({SessionList.group_label}, 'WildType'));
Casp_ind = find(strcmp({SessionList.group_label}, 'IslandKilled'));

% decoding results
Match_M = Acc_Ratio_celltype{celltype_i_temp};
Control = Match_M(Cont_ind,:);
Casp    = Match_M(Casp_ind,:);



%% plot results
celltype = celltype_i_temp;

Cont_color = [1 1 1]*0.2;
Casp_color =[255 168 60]/255;

M_Cont = mean(Control, 'omitmissing');
M_Casp = mean(Casp, 'omitmissing');

SEM_Cont = std(Control, 0, 'omitmissing') ./ sqrt(sum(~isnan(Control)));
SEM_Casp = std(Casp, 0, 'omitmissing') ./ sqrt(sum(~isnan(Casp)));

pos = [65   678   810   275];
fig = figure('Position', pos);
% fig = figure;
% subplot(1,2,1)
hold on

% --- scatter ---
for i = 1:length(CellN)
    jitter = 0.5;
    scatter(CellN(i) - jitter*0.3*randn(size(Control(:,i))), Control(:,i), ...
        10, Cont_color, 'filled', 'MarkerFaceAlpha', 0.3)
    scatter(CellN(i) + jitter*0.3*randn(size(Casp(:,i))), Casp(:,i), ...
        10, Casp_color, 'filled', 'MarkerFaceAlpha', 0.3)
end

plot(CellN, M_Cont, 'Color', Cont_color, 'LineWidth', 3)
plot(CellN, M_Casp, 'Color', Casp_color, 'LineWidth', 3)

errorbar(CellN, M_Cont, SEM_Cont, ...
    'Color', Cont_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)

errorbar(CellN, M_Casp, SEM_Casp, ...
    'Color', Casp_color, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.5)

if b == 1
    ylim([0 0.5])
elseif b == 2
    ylim([0 0.8])
end


xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('decoding accuracy [ratio]')
box off
set(gca,'FontSize',10)



%% create fitting model

[nA_ctrl, nN] = size(Control);
[nA_casp, nN2] = size(Casp);

Accuracy = [];
CellN_all = [];

Group = strings(0);
Session = strings(0);
Animal = strings(0);
RecDay = strings(0);

% Control
for a = 1:nA_ctrl
    for n = 1:nN
        val = Control(a,n);
        if ~isnan(val)
            Accuracy(end+1,1) = val;
            CellN_all(end+1,1) = CellN(n);
            Group(end+1,1) = "Control";
            Session(end+1,1) = "C" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a);
        end
    end
end

% Casp
for a = 1:nA_casp
    for n = 1:nN2
        val = Casp(a,n);
        if ~isnan(val)
            Accuracy(end+1,1) = val;
            CellN_all(end+1,1) = CellN(n);
            Group(end+1,1) = "Casp";
            Session(end+1,1) = "K" + a;

            A = {SessionList.animal_name};
            Animal(end+1,1) = A(a + nA_ctrl);
            R = {SessionList.recording_day};
            RecDay(end+1,1) = R(a + nA_ctrl);
        end
    end
end

% remove invalid rows
valid_idx = ...
    ~isnan(Accuracy) & ...
    ~isnan(CellN_all) & ...
    ~isinf(Accuracy) & ...
    ~isinf(CellN_all);

Accuracy = Accuracy(valid_idx);
CellN_all = CellN_all(valid_idx);
Group = Group(valid_idx);
Animal = Animal(valid_idx);

T = table(Accuracy,CellN_all,categorical(Group),categorical(Animal), ...
    categorical(RecDay),categorical(Session), ...
    'VariableNames',{'Accuracy','N','Group','Animal','RecDay','Session'});
T = rmmissing(T);

lme = fitlme(T,'Accuracy ~ N*Group + (N|Session)');
anova(lme,'DFMethod','Satterthwaite')


%% visualization

Nfit = linspace(min(CellN), max(CellN), 200)';

Tpred_Control = table(Nfit, ...
    categorical(repmat("Control",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'N','Group','Session'});

Tpred_Casp = table(Nfit, ...
    categorical(repmat("Casp",length(Nfit),1)), ...
    categorical(ones(length(Nfit),1)), ...
    'VariableNames',{'N','Group','Session'});

yfit_Control = predict(lme,Tpred_Control);
yfit_Casp = predict(lme,Tpred_Casp);

% figure
subplot(1,2,2)
hold on


MarkerSize = 6;
Alpha = 0.2;

% individual points
for i = 1:length(CellN)
    jitter = 0.6;
    scatter(CellN(i) - jitter*randn(size(Control(:,i))), ...
        Control(:,i), MarkerSize, Cont_color, ...
        'filled','MarkerFaceAlpha', Alpha)
    scatter(CellN(i) + jitter*randn(size(Casp(:,i))), ...
        Casp(:,i), MarkerSize, Casp_color, ...
        'filled','MarkerFaceAlpha', Alpha)
end

% model fit
plot(Nfit,yfit_Control,'Color',Cont_color,'LineWidth',2)
plot(Nfit,yfit_Casp,'Color',Casp_color,'LineWidth',2)

xlabel(strcat("Number of ", celltype_str{celltype}') )
ylabel('decoding accuracy [ratio]')

% ylim([0 0.7])
if b == 1
    % ylim([0 0.5])
    ylim([0 0.4])
    yticks([0 0.2 0.4])
elseif b == 2
    ylim([0 0.8])
end
box off
set(gca,'FontSize',10)

xticks([0 30 60])

anova_tbl = anova(lme,'DFMethod','Satterthwaite'); % 型は table
p_Group = anova_tbl.pValue(strcmp(anova_tbl.Term,'Group'));
p_Int = anova_tbl.pValue(strcmp(anova_tbl.Term,'N:Group'));


%%
%% print statistics
aov = anova(lme,'DFMethod','Satterthwaite');

fprintf('===== Decoding Accuracy =====\n')
fprintf('===== Mixed-effects model =====\n')

% split table
T_ctrl = T(T.Group == "Control", :);
T_casp = T(T.Group == "Casp", :);

% Control
n_obs_ctrl = height(T_ctrl);
n_session_ctrl = numel(unique(T_ctrl.Session));
n_animal_ctrl = numel(unique(T_ctrl.Animal));

% Casp
n_obs_casp = height(T_casp);
n_session_casp = numel(unique(T_casp.Session));
n_animal_casp = numel(unique(T_casp.Animal));

% print
fprintf('\n===== Sample size =====\n')
fprintf('Control: n = %d observations, %d sessions, %d animals\n', ...
    n_obs_ctrl, n_session_ctrl, n_animal_ctrl);
fprintf('Casp:    n = %d observations, %d sessions, %d animals\n', ...
    n_obs_casp, n_session_casp, n_animal_casp);


% stats
for i = 1:height(aov)
    term = string(aov.Term(i));
    F = aov.FStat(i);
    df1 = aov.DF1(i);
    df2 = aov.DF2(i);
    p = aov.pValue(i);
    if p < 0.0001
        pstr = "p < 0.0001";
    else
        pstr = sprintf("p = %.4g", p);
    end
    fprintf('%s: %s, F (%.3f, %.3f) = %.3f\n', ...
        term, pstr, df1, df2, F);
end

end


