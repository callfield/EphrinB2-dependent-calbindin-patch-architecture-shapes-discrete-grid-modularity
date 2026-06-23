clear
close all

%% output folder

OutDir = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\stats_260206\StatsResults_260524\J6C7_hist_noP1';
[status,msg,msgID] = mkdir(fullfile(OutDir)); %#ok<ASGLU>

OutDir_panel = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\stats_260206\StatsResults_260524\J6C7_hist_noP1\ForFigPanel';
[status,msg,msgID] = mkdir(fullfile(OutDir_panel)); %#ok<ASGLU>


[status,msg,msgID] = mkdir(fullfile(OutDir_panel, 'Field')); %#ok<ASGLU>
[status,msg,msgID] = mkdir(fullfile(OutDir, 'Field')); %#ok<ASGLU>


%% load data

% load all mice data
mice_str = {'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'n12', 'n14', ...
    'p1', 'p7','p8', 'p11', 'p13', 'p14','p15', 'p16'};
mice_str_Folder = {'N5', 'N6', 'N7', 'N8', 'N9', 'N10', 'N12', 'N14', ...
    'P1', 'P7','P8', 'P11', 'P13', 'P14','P15', 'P16'};


Nanimal = 16;
Ndays = 7;

N_NegAnimal = 8; % Number of Control mice
N_PosAnimal = 8; % Number of Casp mice

LoadData = cell(Nanimal, Ndays);

for Day = 1:7
    for mice_ind = 1:Nanimal
        if (mice_ind >= 13 || (mice_ind >=7 && mice_ind<=8) )
            Dir = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\CleanedAnalysisData\CleanedAnalysisData_AdpBin_ThrPeakRatio';
        else
            Dir = 'H:\experiments H drive\250123 Ca imaging\Matlab Codes\Place Field analysis\statistics\CleanedAnalysisData\CleanedAnalysisData_AdpBin_ThrPeakRatio';
        end

        s_savemat = strcat('CleanedAnalysisData_', mice_str{mice_ind}, '_d', num2str(Day), '.mat');
        ss = fullfile(Dir, s_savemat);

        fprintf(strcat("load ", s_savemat, "\n") );
        if exist(ss, 'file') == 0
            fprintf(strcat("File not found: ", 'Mouse', mice_str_Folder{mice_ind}, " Day%d\n"), Day);
            continue;
        end
        LoadData{mice_ind, Day} = load(ss);

    end
end




%% summarize data

close all
Boot = cell(Nanimal, Ndays);

thr_p = 0.05;
thr_delta = 0.147;

CellListT_all_animals = cell(Nanimal, Ndays);
FieldStats_table_animal = cell(Nanimal, Ndays);
Cell_Cetegory_Day = cell(1, Ndays);
CellList_category_all = cell(Nanimal, Ndays);
shuf_CellList_category_all = cell(Nanimal, Ndays);
RunDist = nan(Nanimal, Ndays);
RunDist_2 = nan(Nanimal, Ndays);
RunDurationRatio = nan(Nanimal, Ndays);
RunDurationRatio_2 = nan(Nanimal, Ndays);

Velocity_2 = nan(Nanimal, Ndays);

CellListT_all_animals_shuf = cell(Nanimal, Ndays);


% Remove if running coverage does not satisfy threshold
RunCheck = load("H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\stats_260206\OccRate\Occresults_raw_260208.mat");
RunCoverageCheck = (RunCheck.Jaccard > 0.6) & (RunCheck.Coverage_All > 0.7);

for Day = 1:1:7

    Cell_Cetegory = nan(Nanimal,4);

    for mice_ind = 1:Nanimal
        % Basics

        % Cell_Cetegory(mice_ind,:) = [NaN NaN NaN NaN];
        if isempty(LoadData{mice_ind, Day} )
            continue;
        end
        RunDist(mice_ind, Day) = sum(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed, "omitmissing") * (1/LoadData{mice_ind, Day}.BasicAnalysis.VideoSampling) / 100; %ran distance [m]

        RunDurationRatio(mice_ind, Day) = length( find(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed >2 ) ) / length(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed  ) ;

        % Remove if running coverage does not satisfy RunCoverageCheck = (RunCheck.Jaccard > 0.6) & (RunCheck.Coverage_All > 0.7);

        % N14 is removed because of circling behavior
        % P1 and P13 is removed because of failure of apoosis inducing at ipsilateral MEC with recording site (right MEC)
        RunCoverageCheck(8,:) = 0; %N14
        RunCoverageCheck(13,:) = 0; %P13
        RunCoverageCheck(9,:) = 0; %P1

        if RunCoverageCheck(mice_ind,Day) == 0
            continue;
        end


        RunDist_2(mice_ind, Day) = sum(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed, "omitmissing") * (1/LoadData{mice_ind, Day}.BasicAnalysis.VideoSampling) / 100; %ran distance [m]
        RunDurationRatio_2(mice_ind, Day) = length( find(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed >2 ) ) / length(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed  ) ;

        Velocity_2(mice_ind, Day) = mean(LoadData{mice_ind, Day}.BasicAnalysis.Detections.AnimalSpeed, 'omitmissing');

        Boot{mice_ind, Day} = LoadData{mice_ind, Day}.BOOT;

        [CellFieldList_table, CellFieldList_shuf_table, CellSpeedModList_table] = func_GenerateCellList(LoadData, Boot, mice_ind, Day);

        CellList_table = join(CellFieldList_table, CellSpeedModList_table);
        CellListT_all_animals{mice_ind, Day} = CellList_table;
        CellListT_all_animals_shuf{mice_ind, Day} = CellFieldList_shuf_table;

        %% Field

        FieldStats_table = func_GenerateFieldList(LoadData, mice_ind, Day, CellList_table);
        FieldStats_table_animal{mice_ind, Day} = FieldStats_table;

        %% Define cell categories
        Ind_uncategorized = CellList_table.Rel_boot_p >= thr_p | abs(CellList_table.Rel_boot_delta) < thr_delta | isnan(CellList_table.Rel_boot_p);
        Cell_uncategorized = CellList_table(Ind_uncategorized,:);

        Ind_unstableRatemap =  CellList_table.Rel_boot_p < thr_p & CellList_table.Rel_boot_delta <= -thr_delta;
        Cell_unstableRatemap = CellList_table(Ind_unstableRatemap,:);

        Ind_reliableRatemap =  CellList_table.Rel_boot_p < thr_p & CellList_table.Rel_boot_delta >= thr_delta;
        Cell_reliableRatemap = CellList_table(Ind_reliableRatemap,:);

        Ind_reliableRatemap_PlaceCell =  CellList_table.Rel_boot_p < thr_p & CellList_table.Rel_boot_delta >= thr_delta & CellList_table.SpInfo_p < 0.05;
        Cell_reliableRatemap_PlaceCell = CellList_table(Ind_reliableRatemap_PlaceCell,:);

        Ind_reliableRatemap_NotPlaceCell =  CellList_table.Rel_boot_p < thr_p & CellList_table.Rel_boot_delta >= thr_delta & CellList_table.SpInfo_p >= 0.05;
        Cell_reliableRatemap_NotPlaceCell = CellList_table(Ind_reliableRatemap_NotPlaceCell,:);

        Cell_Cetegory(mice_ind,:) = [height(Cell_uncategorized), height(Cell_unstableRatemap), height(Cell_reliableRatemap_NotPlaceCell), height(Cell_reliableRatemap_PlaceCell)];


        CellList_category.Cell_uncategorized = Cell_uncategorized;
        CellList_category.Cell_unstableRatemap = Cell_unstableRatemap;
        CellList_category.Cell_reliableRatemap = Cell_reliableRatemap;
        CellList_category.Cell_reliableRatemap_PlaceCell = Cell_reliableRatemap_PlaceCell;
        CellList_category.Cell_reliableRatemap_NotPlaceCell = Cell_reliableRatemap_NotPlaceCell;

        CellList_category_all{mice_ind, Day} = CellList_category;



        %% for shuf
        shuf_Cell_uncategorized                  = CellFieldList_shuf_table(Ind_uncategorized,:);
        shuf_Cell_unstableRatemap                = CellFieldList_shuf_table(Ind_unstableRatemap,:);
        shuf_Cell_reliableRatemap                = CellFieldList_shuf_table(Ind_reliableRatemap,:);
        shuf_Cell_reliableRatemap_PlaceCell      = CellFieldList_shuf_table(Ind_reliableRatemap_PlaceCell,:);
        shuf_Cell_reliableRatemap_NotPlaceCell   = CellFieldList_shuf_table(Ind_reliableRatemap_NotPlaceCell,:);

        shuf_Cell_Cetegory(mice_ind,:) = [height(shuf_Cell_uncategorized), height(shuf_Cell_unstableRatemap), height(shuf_Cell_reliableRatemap_NotPlaceCell), height(shuf_Cell_reliableRatemap_PlaceCell)];


        shuf_CellList_category.Cell_uncategorized = shuf_Cell_uncategorized;
        shuf_CellList_category.Cell_unstableRatemap = shuf_Cell_unstableRatemap;
        shuf_CellList_category.Cell_reliableRatemap = shuf_Cell_reliableRatemap;
        shuf_CellList_category.Cell_reliableRatemap_PlaceCell = shuf_Cell_reliableRatemap_PlaceCell;
        shuf_CellList_category.Cell_reliableRatemap_NotPlaceCell = shuf_Cell_reliableRatemap_NotPlaceCell;

        shuf_CellList_category_all{mice_ind, Day} = shuf_CellList_category;


    end

    varNames2 = {'uncategorized', 'unstable', 'non place cell', 'place cell'};
    Cell_Cetegory_table = array2table(Cell_Cetegory, 'VariableNames',varNames2);
    Cell_Cetegory_Day{Day} = Cell_Cetegory_table;

end



%% Running behavior

close all

RunDist = RunDist_2;

nDay = 7;
fig = figure('Position', [-1584         550        700         220]);
sgst = strcat("Running distance");
sgtitle(sgst)
subplot(1,7,[1 3])
hold on
A = RunDist(1:N_NegAnimal,:);
B = RunDist(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

errorbar(1:nDay, mA, semA, '-o', 'Color',[0 0 1], 'LineWidth',1.5, ...
    'MarkerFaceColor',[0 0 1], 'CapSize',8);
errorbar(1:nDay, mB, semB, '-o', 'Color',[1 0 0], 'LineWidth',1.5, ...
    'MarkerFaceColor',[1 0 0], 'CapSize',8);

xlim([0.5 nDay+0.5]);
xticks(1:nDay);
xlabel('Day');
ylabel('Running distance [m]');
ylim([0 250])
legend({'Negative','Positive'}, 'Location','best');

yL = ylim;
for i = 1:nDay
    [h,p,ci,stats] = ttest2(A(:, i), B(:, i));

    yPos = yL(2) - 0.05*range(yL);
    if p < 0.001
        text(i, yPos, '***', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    elseif p <= 0.01
        text(i, yPos, '**', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    elseif p < 0.05
        text(i, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    end
    % disp(p)
end
title('Day time lapse')
set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])


%
for i = 1:4
    subplot(1,7,3+i)
    hold on

    if i == 1
        A = mean(RunDist(1:N_NegAnimal,:), 2,"omitmissing");
        B = mean(RunDist(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:), 2,"omitmissing");
        title('Day all')
    elseif i == 2
        A = mean(RunDist(1:N_NegAnimal, 1:3), 2,"omitmissing");
        B = mean(RunDist(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,1:3), 2,"omitmissing");
        title('Day 1-3')
    elseif i == 3
        A = mean(RunDist(1:N_NegAnimal, 4:7), 2,"omitmissing");
        B = mean(RunDist(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,4:7), 2,"omitmissing");
        title('Day 4-7')

    elseif i == 4 % first three successful trials

        A = RunDist(1:N_NegAnimal,:);
        B = RunDist(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);

        [A, ind_NegNonNaN] = first3nonNaN(A);
        [B, ind_PosNonNaN] = first3nonNaN(B);

        title('3 successful trials')
    end
    A = A(:);
    B = B(:);

    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
    semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

    errorbar(1, mA, semA, 'o', 'Color',[0 0 1], 'LineWidth',1.5, ...
        'MarkerFaceColor',[0 0 1], 'CapSize',8);
    errorbar(2, mB, semB, 'o', 'Color',[1 0 0], 'LineWidth',1.5, ...
        'MarkerFaceColor',[1 0 0], 'CapSize',8);

    xlim([0.5 2.5]);
    xticks(1:2);
    xticklabels({'Neg', 'Pos'})
    ylim([0 250])

    [h,p,ci,stats] = ttest2(A, B);
    yL = ylim;
    yPos = yL(2) - 0.05*range(yL);
    if p < 0.05
        text(1.5, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    else
        text(1.5, yPos, 'n.s.', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    end
    hold off;

    set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])
end

s = strcat("Running distance.jpg");
exportgraphics(fig, fullfile(OutDir, s));

%%

A = RunDist_2(1:N_NegAnimal,:);
B = RunDist_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);

[A, ind_NegNonNaN] = first3nonNaN(A);
[B, ind_PosNonNaN] = first3nonNaN(B);
title('3 successful trials')

figure
hold on

nAnimalA = size(A,1);
nTrialA  = size(A,2);

nAnimalB = size(B,1);
nTrialB  = size(B,2);

colorsA = lines(nAnimalA);
colorsB = lines(nAnimalB);

jitterWidth = 0.15;

% ----- jitter -----
jitA = (rand(nAnimalA,1)-0.5)*jitterWidth;
jitB = (rand(nAnimalB,1)-0.5)*jitterWidth;

% ----- Neg animals -----
for a = 1:nAnimalA

    xAnimal = 1 + jitA(a);

    for t = 1:nTrialA

        if ~isnan(A(a,t))

            scatter(xAnimal, A(a,t), 25, ...
                'MarkerEdgeColor', colorsA(a,:), ...
                'MarkerFaceColor', colorsA(a,:), ...
                'MarkerFaceAlpha', 0.5, ...
                'MarkerEdgeAlpha', 0.5);
        end
    end
end


% ----- Pos animals -----
for a = 1:nAnimalB

    xAnimal = 2 + jitB(a);

    for t = 1:nTrialB

        if ~isnan(B(a,t))

            scatter(xAnimal, B(a,t), 25, ...
                'MarkerEdgeColor', colorsB(a,:), ...
                'MarkerFaceColor', colorsB(a,:), ...
                'MarkerFaceAlpha', 0.5, ...
                'MarkerEdgeAlpha', 0.5);
        end
    end
end


meanA = mean(A,2,"omitmissing");
meanB = mean(B,2,"omitmissing");

semAnimalA = std(A,0,2,"omitmissing") ./ sqrt(sum(~isnan(A),2));
semAnimalB = std(B,0,2,"omitmissing") ./ sqrt(sum(~isnan(B),2));

xA = 1 + jitA;
xB = 2 + jitB;

errorbar(xA, meanA, semAnimalA, 'ko', ...
    'MarkerFaceColor','w', ...
    'CapSize',4, ...
    'LineWidth',1);

errorbar(xB, meanB, semAnimalB, 'ko', ...
    'MarkerFaceColor','w', ...
    'CapSize',4, ...
    'LineWidth',1);


% average plot
A = A(:);
B = B(:);

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

errorbar(1, mA, semA, 'o', 'Color',[0 0 1], 'LineWidth',1.5, ...
    'MarkerFaceColor',[0 0 1], 'CapSize',8);
errorbar(2, mB, semB, 'o', 'Color',[1 0 0], 'LineWidth',1.5, ...
    'MarkerFaceColor',[1 0 0], 'CapSize',8);

xlim([0.5 2.5]);
xticks(1:2);
xticklabels({'Neg', 'Pos'})
ylim([0 250])

[h,p,ci,stats] = ttest2(A, B);
yL = ylim;
yPos = yL(2) - 0.05*range(yL);
if p < 0.05
    text(1.5, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
else
    text(1.5, yPos, 'n.s.', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
end
hold off;

set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])


%% Running Duration ratio

close all

nDay = 7;


fig = figure('Position', [-1584         550        700         220]);
sgst = strcat("Running duration");
sgtitle(sgst)
subplot(1,7,[1 3])
hold on

A = RunDurationRatio_2(1:N_NegAnimal,:) * 100;
B = RunDurationRatio_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:) * 100;

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

errorbar(1:nDay, mA, semA, '-o', 'Color',[0 0 1], 'LineWidth',1.5, ...
    'MarkerFaceColor',[0 0 1], 'CapSize',8);
errorbar(1:nDay, mB, semB, '-o', 'Color',[1 0 0], 'LineWidth',1.5, ...
    'MarkerFaceColor',[1 0 0], 'CapSize',8);

xlim([0.5 nDay+0.5]);
xticks(1:nDay);
xlabel('Day');
ylabel('Running distance [m]');
ylim([0 50])
legend({'Negative','Positive'}, 'Location','best');

yL = ylim;
for i = 1:nDay
    [h,p,ci,stats] = ttest2(A(:, i), B(:, i));

    yPos = yL(2) - 0.05*range(yL);
    if p < 0.001
        text(i, yPos, '***', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    elseif p <= 0.01
        text(i, yPos, '**', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    elseif p < 0.05
        text(i, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    end
end

title('Day time lapse')
set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])


%
for i = 1:4
    subplot(1,7,3+i)
    hold on
    if i == 1
        A = mean(RunDurationRatio_2(1:N_NegAnimal,:), 2,"omitmissing");
        B = mean(RunDurationRatio_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:), 2,"omitmissing");
        title('Day all')
    elseif i == 2
        A = mean(RunDurationRatio_2(1:N_NegAnimal, 1:3), 2,"omitmissing");
        B = mean(RunDurationRatio_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal, 1:3), 2,"omitmissing");
        title('Day 1-3')
    elseif i == 3
        A = mean(RunDurationRatio_2(1:N_NegAnimal, 1:3), 2,"omitmissing");
        B = mean(RunDurationRatio_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal, 1:3), 2,"omitmissing");
        title('Day 4-7')
    elseif i == 4 % first three successful trials
        A = RunDurationRatio_2(1:N_NegAnimal,:);
        B = RunDurationRatio_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);
        A = first3nonNaN(A);
        B = first3nonNaN(B);

        title('3 successful trials')
    end
    A = A(:) *100;
    B = B(:) *100;

    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
    semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

    errorbar(1, mA, semA, 'o', 'Color',[0 0 1], 'LineWidth',1.5, ...
        'MarkerFaceColor',[0 0 1], 'CapSize',8);
    errorbar(2, mB, semB, 'o', 'Color',[1 0 0], 'LineWidth',1.5, ...
        'MarkerFaceColor',[1 0 0], 'CapSize',8);

    xlim([0.5 2.5]);
    xticks(1:2);
    xticklabels({'Neg', 'Pos'})
    % ylim([0 250])
    ylim([0 50])

    [h,p,ci,stats] = ttest2(A, B);
    yL = ylim;
    yPos = yL(2) - 0.05*range(yL);
    if p < 0.05
        text(1.5, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    else
        text(1.5, yPos, 'n.s.', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    end
    hold off;

    set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])
end


%% For Fig
Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

close all

fig = figure;
hold on
set(fig, 'Position', [651   726   120   260])

A = RunDurationRatio_2(1:N_NegAnimal,:);
B = RunDurationRatio_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);
A = first3nonNaN(A);
B = first3nonNaN(B);
A = mean(A,2,"omitmissing");
B = mean(B,2,"omitmissing");

title('Day 1-3')

A = A(:) *100;
B = B(:) *100;

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

hold on
bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
ax = gca;
ax.Layer = 'top';

errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',8);
errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',8);


% ---- smart jitter ----
width = 0.08;
nbins = 20;

edgesA = linspace(min(A), max(A) + eps(max(A)), nbins+1);
edgesB = linspace(min(B), max(B) + eps(max(B)), nbins+1);
xA = nan(size(A));
for j = 1:nbins
    idx = A >= edgesA(j) & A < edgesA(j+1);
    n = sum(idx);
    if n == 1
        xA(idx) = 1;
    elseif n > 1
        offsets = linspace(-width, width, n);
        xA(idx) = 1 + offsets;
    end
end

xB = nan(size(B));
for j = 1:nbins
    idx = B >= edgesB(j) & B < edgesB(j+1);
    n = sum(idx);
    if n == 1
        xB(idx) = 2;
    elseif n > 1
        offsets = linspace(-width, width, n);
        xB(idx) = 2 + offsets;
    end
end
plot(xA, A, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth',0.5);
plot(xB, B, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth',0.5);


xticks(1:2);
xticklabels({'Neg', 'Pos'})


[h,p,ci,stats] = ttest2(A, B);
yL = ylim;
yPos = yL(2) - 0.05*range(yL);
hold off;

set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])


ylim([0 50])
yticks([0 25 50])
yticklabels([])
xticklabels([])
title('')
xlim([0.4 2.6])


s = strcat("Running duration.emf");
exportgraphics(fig, fullfile(OutDir_panel, s));

sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('A: %.2f ± %.2f\n', mA, sdA);
fprintf('B: %.2f ± %.2f\n', mB, sdB);

%% mean speed

close all

nDay = 7;

fig = figure('Position', [-1584         550        700         220]);
sgst = strcat("Running duration");
sgtitle(sgst)
subplot(1,7,[1 3])
hold on
A = Velocity_2(1:N_NegAnimal,:);
B = Velocity_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);


mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

errorbar(1:nDay, mA, semA, '-o', 'Color',[0 0 1], 'LineWidth',1.5, ...
    'MarkerFaceColor',[0 0 1], 'CapSize',8);
errorbar(1:nDay, mB, semB, '-o', 'Color',[1 0 0], 'LineWidth',1.5, ...
    'MarkerFaceColor',[1 0 0], 'CapSize',8);

xlim([0.5 nDay+0.5]);
xticks(1:nDay);
xlabel('Day');
ylabel('Mean Velocity [cm/s]');
legend({'Negative','Positive'}, 'Location','best');

yL = ylim;
for i = 1:nDay
    [h,p,ci,stats] = ttest2(A(:, i), B(:, i));

    yPos = yL(2) - 0.05*range(yL);
    if p < 0.001
        text(i, yPos, '***', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    elseif p <= 0.01
        text(i, yPos, '**', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    elseif p < 0.05
        text(i, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    end
end
title('Day time lapse')
set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])


%
for i = 1:4
    subplot(1,7,3+i)
    hold on

    if i == 1
        A = mean(Velocity_2(1:N_NegAnimal,:), 2,"omitmissing");
        B = mean(Velocity_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:), 2,"omitmissing");

        title('Day all')
    elseif i == 2
        A = mean(Velocity_2(1:N_NegAnimal, 1:3), 2,"omitmissing");
        B = mean(Velocity_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal, 1:3), 2,"omitmissing");

        title('Day 1-3')
    elseif i == 3
        A = mean(Velocity_2(1:N_NegAnimal, 1:3), 2,"omitmissing");
        B = mean(Velocity_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal, 1:3), 2,"omitmissing");

        title('Day 4-7')
    elseif i == 4
        A = Velocity_2(1:N_NegAnimal,:) ;
        B = Velocity_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);

        A = first3nonNaN(A);
        B = first3nonNaN(B);
        A = mean(A,2,"omitmissing");
        B = mean(B,2,"omitmissing");

        title('3 successful trials')
    end
    A = A(:) ;
    B = B(:) ;

    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
    semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

    errorbar(1, mA, semA, 'o', 'Color',[0 0 1], 'LineWidth',1.5, ...
        'MarkerFaceColor',[0 0 1], 'CapSize',8);
    errorbar(2, mB, semB, 'o', 'Color',[1 0 0], 'LineWidth',1.5, ...
        'MarkerFaceColor',[1 0 0], 'CapSize',8);

    xlim([0.5 2.5]);
    xticks(1:2);
    xticklabels({'Neg', 'Pos'})

    [h,p,ci,stats] = ttest2(A, B);
    yL = ylim;
    yPos = yL(2) - 0.05*range(yL);
    if p < 0.05
        text(1.5, yPos, '*', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    else
        text(1.5, yPos, 'n.s.', 'FontSize', 12, 'HorizontalAlignment','center', 'Color','k');
    end
    hold off;

    set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])
end

%% Speed bar plot
% For Fig
Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

close all
% For Fig
fig = figure;
hold on
set(fig, 'Position', [651   726   120   260])


A = Velocity_2(1:N_NegAnimal,:) ;
B = Velocity_2(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:);

A = first3nonNaN(A);
B = first3nonNaN(B);
A = mean(A,2,"omitmissing");
B = mean(B,2,"omitmissing");

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

hold on
bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
ax = gca;
ax.Layer = 'top';

errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',8);
errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',8);

% ---- smart jitter ----
width = 0.08;
nbins = 20;
edgesA = linspace(min(A), max(A) + eps(max(A)), nbins+1);
edgesB = linspace(min(B), max(B) + eps(max(B)), nbins+1);
xA = nan(size(A));
for j = 1:nbins
    idx = A >= edgesA(j) & A < edgesA(j+1);
    n = sum(idx);
    if n == 1
        xA(idx) = 1;
    elseif n > 1
        offsets = linspace(-width, width, n);
        xA(idx) = 1 + offsets;
    end
end

xB = nan(size(B));
for j = 1:nbins
    idx = B >= edgesB(j) & B < edgesB(j+1);
    n = sum(idx);
    if n == 1
        xB(idx) = 2;
    elseif n > 1
        offsets = linspace(-width, width, n);
        xB(idx) = 2 + offsets;
    end
end

plot(1, A, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth',0.5);
plot(2, B, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'LineWidth',0.5);

xticks(1:2);
xticklabels({'Neg', 'Pos'})


[h,p,ci,stats] = ttest2(A, B);
yL = ylim;
yPos = yL(2) - 0.05*range(yL);

hold off;

set(gca, 'Position', get(gca,'Position') + [0 0.05 0 -0.15])

ylim([0 4])
yticks([0 2 4])
yticklabels([])
xticklabels([])
title('')
xlim([0.4 2.6])

s = strcat("Velocity.emf");
exportgraphics(fig, fullfile(OutDir_panel, s));

sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('A: %.2f ± %.2f\n', mA, sdA);
fprintf('B: %.2f ± %.2f\n', mB, sdB);


%% cell categories

close all

pval = nan(10,1);
p_fisher = nan(10,1);
Cell_Cetegory_acrossDay = nan(10,7);
for Day = 11 %1:10

    if  Day == 9
        Cell_Cetegory_temp = ( cellfun(@(x) table2array(x), Cell_Cetegory_Day, 'UniformOutput', false) )';
        sumMat = sum(cat(3, Cell_Cetegory_temp{1:3}), 3, "omitmissing");
    end

    %
    if Day == 11 % successful three trials
        Cell_Cetegory_temp = ( cellfun(@(x) table2array(x), Cell_Cetegory_Day, 'UniformOutput', false) )';
        nAnimal = 16;
        nDay = 7;
        First3Days = NaN(nAnimal, 4, 3);
        First3Days_cell = cell(3,1);

        for a = 1:nAnimal
            validDays = [];
            % ---- find successful recording day ----
            for d = 1:nDay
                data = Cell_Cetegory_temp{d};   % 16×4

                if ~all(isnan(data(a,:)))       % find non-nan animals
                    validDays(end+1) = d;
                end
                if numel(validDays) == 3
                    break
                end
            end

            % ---- data extraction ----
            for k = 1:numel(validDays)
                d = validDays(k);
                First3Days(a,:,k) = Cell_Cetegory_temp{d}(a,:);
                %
            end
        end

        for k = 1:3
            First3Days_cell{k} = First3Days(:,:,k);
        end

        sumMat = zeros(nAnimal,4);
        for k = 1:3
            tmp = First3Days_cell{k};
            tmp(isnan(tmp)) = 0;
            sumMat = sumMat + tmp;
        end
    end

    Cell_Cetegory = [sumMat(:,1) + sumMat(:,2), sumMat(:,3), sumMat(:,4)];

    Color = nan(4,3);
    Color(1, 1:3) = [1, 1, 1]*0.8;
    Color(2, 1:3) = [0.7, 0.4, 1];
    Color(3, 1:3) = [1, 0.0, 0];

    fig1 = figure('Position', [-1584         404        1486         221]);
    subplot(151)
    hold on
    h_bar = bar(1:Nanimal, Cell_Cetegory, 'stacked');
    h_bar(1).FaceColor = Color(1, :); h_bar(2).FaceColor = Color(2, :); h_bar(3).FaceColor = Color(3, :); %h_bar(4).FaceColor = Color(4, :);

    ax = gca;
    ax.XTick = 1:Nanimal; ax.XTickLabel = mice_str_Folder;
    ylabel('Cell Counts')
    hold on
    varNames2 = {'uncategorized', 'unstable', 'non place cell', 'place cell'};

    subplot(152)
    hold on
    h_bar = bar(1:Nanimal, Cell_Cetegory ./ sum(Cell_Cetegory, 2) *100, 'stacked');
    h_bar(1).FaceColor = Color(1, :); h_bar(2).FaceColor = Color(2, :); h_bar(3).FaceColor = Color(3, :); %h_bar(4).FaceColor = Color(4, :);
    ax = gca;
    ax.XTick = 1:Nanimal; ax.XTickLabel = mice_str_Folder;
    ylabel('Cell Count Percentage')
    ylim([0 100])

    subplot(153)
    hold on
    A = Cell_Cetegory(:,2:3) ./ sum(Cell_Cetegory(:,2:3), 2) *100;
    h_bar = bar(1:Nanimal, A, 'stacked');
    h_bar(1).FaceColor = Color(2, :); h_bar(2).FaceColor = Color(3, :);
    ax = gca;
    ax.XTick = 1:Nanimal; ax.XTickLabel = mice_str_Folder;
    ylabel('Place cell Percentage')
    ylim([0 100])

    if Day<=7
        Cell_Cetegory_acrossDay(:,Day) = A(:,2);  % percentage of place cell
    end

    subplot(154)
    hold on
    A = sum(Cell_Cetegory(1:N_NegAnimal,:), 'omitmissing');
    B = sum(Cell_Cetegory(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:), 'omitmissing');
    Arate = A / sum(A) *100;
    Brate = B / sum(B) *100;
    h_bar = bar(1:2, [Arate; Brate], 'stacked');
    h_bar(1).FaceColor = Color(1, :); h_bar(2).FaceColor = Color(2, :); h_bar(3).FaceColor = Color(3, :); %h_bar(4).FaceColor = Color(4, :);
    ax = gca;
    ax.XTick = 1:2; ax.XTickLabel = {'N', 'P'};
    ylabel('Cell Count Percentage')
    ylim([0 100])

    subplot(155)
    hold on
    A = sum(Cell_Cetegory(1:N_NegAnimal,2:3), 'omitmissing');
    B = sum(Cell_Cetegory(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,2:3), 'omitmissing');
    Arate = A / sum(A) *100;
    Brate = B / sum(B) *100;

    h_bar = bar(1:2, [Arate; Brate], 'stacked');
    h_bar(1).FaceColor = Color(2, :); h_bar(2).FaceColor = Color(3, :);
    ax = gca;
    ax.XTick = 1:2; ax.XTickLabel = {'N', 'P'};
    ylim([0 100])

    % chi2 test
    tbl_cross = [A; B];   % cross table
    rowSum = sum(tbl_cross, 2);
    colSum = sum(tbl_cross, 1);
    n = sum(tbl_cross, 'all');
    expected = rowSum * colSum / n;
    chi2stat = sum((tbl_cross - expected).^2 ./ expected, 'all');
    df = (size(tbl_cross,1)-1) * (size(tbl_cross,2)-1);
    pval(Day) = 1 - chi2cdf(chi2stat, df);
    fprintf('Chi2 = %.4f, df = %d, p = %.6f\n', chi2stat, df, pval(Day));

    [h, p_fisher(Day)] = fishertest(tbl_cross);
    fprintf('Fisher exact p = %.6f\n', p_fisher(Day));

    ylabel(sprintf('Place cell Percentage, p = %.4f', pval(Day)) )

    if Day<=7
        sgst = strcat("Day", num2str(Day));
    elseif Day == 8
        sgst = strcat("Day All");
    elseif Day == 9
        sgst = strcat("Day 1-3");
    elseif Day == 10
        sgst = strcat("Day 4-7");
    elseif Day == 11
        sgst = strcat("successful 3 days");
    end
    sgtitle(sgst)

    %%
    close all
    fig = figure;

    set(fig, 'Position', [1099         723         141         155])

    hold on
    A = sum(Cell_Cetegory(1:N_NegAnimal,:), 'omitmissing');
    B = sum(Cell_Cetegory(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:), 'omitmissing');
    Arate = A / sum(A) *100;
    Brate = B / sum(B) *100;
    h_bar = bar(1:2, [Arate; Brate], 'stacked');
    h_bar(1).FaceColor = Color(1, :); h_bar(2).FaceColor = Color(2, :); h_bar(3).FaceColor = Color(3, :);
    for i = 1:3
        h_bar(i).BarWidth = 0.5;
    end
    ax = gca;
    ax.XTick = 1:2; ax.XTickLabel = {'N', 'P'};

    ylim([0 100])
    yticks([0 50 100])
    yticklabels([])
    xticklabels([])
    title('')
    xlim([0.2 2.8])

    % save figs
    s = strcat("CellCounts.emf");
    exportgraphics(fig, fullfile(OutDir_panel, s));


    %% Place cell percentage_WithinReliableCell
    close all
    fig = figure;
    % get(gcf,'Position')
    set(fig, 'Position', [1120         708         120         170])

    hold on
    AA = Cell_Cetegory(:,2:3) ./ sum(Cell_Cetegory(:,2:3), 2) *100;

    A = AA(1:N_NegAnimal,2);
    B = AA(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,2);

    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
    semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

    hold on
    bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
    bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
    ax = gca;
    ax.Layer = 'top';

    errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);
    errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);

    % ---- smart jitter ----
    width = 0.08;
    nbins = 20;
    edgesA = linspace(min(A), max(A), nbins+1);
    xA = nan(size(A));
    for j = 1:nbins
        if j < nbins
            idx = A >= edgesA(j) & A < edgesA(j+1);
        else
            idx = A >= edgesA(j) & A <= edgesA(j+1);
        end
        n = sum(idx);
        if n <= 1
            xA(idx) = 1;
        elseif n > 1
            offsets = linspace(-width, width, n);
            xA(idx) = 1 + offsets;
        end
    end
    edgesB = linspace(min(B), max(B), nbins+1);
    xB = nan(size(B));
    for j = 1:nbins
        if j < nbins
            idx = B >= edgesB(j) & B < edgesB(j+1);
        else
            idx = B >= edgesB(j) & B <= edgesB(j+1);
        end
        n = sum(idx);
        if n <= 1
            xB(idx) = 2;
        elseif n > 1
            offsets = linspace(-width, width, n);
            xB(idx) = 2 + offsets;
        end
    end
    plot(xA, A, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'LineWidth',0.5);
    plot(xB, B, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'LineWidth',0.5);

    xticks(1:2);
    xticklabels({'Neg', 'Pos'})

    [h,p,ci,stats] = ttest2(A, B);
    yL = ylim;
    yPos = yL(2) - 0.05*range(yL);


    ylim([0 100])
    yticks([0 50 100])

    yticklabels([])
    xticklabels([])
    title('')
    xlim([0.2 2.8])

    % save figs
    s = strcat("Place cell percentage_WithinReliableCell.emf");
    exportgraphics(fig, fullfile(OutDir_panel, s));


    disp(stats)

    sdA = std(A,'omitmissing');
    sdB = std(B,'omitmissing');
    fprintf('Control: %.2f ± %.2f\n', mA, sdA);
    fprintf('Casp: %.2f ± %.2f\n', mB, sdB);


    NCont = sum(sum(Cell_Cetegory(1:N_NegAnimal,2:3), 2));
    NCasp = sum(sum(Cell_Cetegory(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,2:3), 2));

end



%% place cell % in all cell
close all
fig = figure;
% get(gcf,'Position')
set(fig, 'Position', [1120         708         120         170])

hold on
AA = Cell_Cetegory ./ sum(Cell_Cetegory, 2) *100;
A = AA(1:N_NegAnimal,2);
B = AA(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,2);

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

hold on
bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
ax = gca;
ax.Layer = 'top';

errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);
errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);

% ---- smart jitter ----
width = 0.08;
nbins = 20;

edgesA = linspace(min(A), max(A) + eps(max(A)), nbins+1);
edgesB = linspace(min(B), max(B) + eps(max(B)), nbins+1);
xA = nan(size(A));
for j = 1:nbins
    if j < nbins
        idx = A >= edgesA(j) & A < edgesA(j+1);
    else
        idx = A >= edgesA(j) & A <= edgesA(j+1);
    end
    n = sum(idx);
    if n <= 1
        xA(idx) = 1;
    elseif n > 1
        offsets = linspace(-width, width, n);
        xA(idx) = 1 + offsets;
    end
end

xB = nan(size(B));
for j = 1:nbins
    if j < nbins
        idx = B >= edgesB(j) & B < edgesB(j+1);
    else
        idx = B >= edgesB(j) & B <= edgesB(j+1);
    end
    n = sum(idx);
    if n <= 1
        xB(idx) = 2;
    elseif n > 1
        offsets = linspace(-width, width, n);
        xB(idx) = 2 + offsets;
    end
end
plot(xA, A, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'LineWidth',0.5);
plot(xB, B, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'w', 'LineWidth',0.5);

xticks(1:2);
xticklabels({'Neg', 'Pos'})

[h,p,ci,stats] = ttest2(A, B);
yL = ylim;
yPos = yL(2) - 0.05*range(yL);

ylim([0 50])
yticks([0 25 50])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])

% save figs
s = strcat("Place cell percentage_InAll.emf");
exportgraphics(fig, fullfile(OutDir_panel, s));

disp(mA)
disp(mB)
disp(semA)
disp(semB)

sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Control: %.2f ± %.2f\n', mA, sdA);
fprintf('Casp: %.2f ± %.2f\n', mB, sdB);

%% CellCategory_AcrossDay (All cells)

close all
fig = figure;
% get(gcf,'Position')
set(fig, 'Position', [965   240   498   345]);

hold on
for mice_ind = 1:10
    p = plot(Cell_Cetegory_acrossDay(mice_ind,:), '-o');
    if mice_ind<=6
        p.Color = 'b';
    else
        p.Color = 'r';
    end
    p.Color(4) = 0.2;
end
ylim([0 100])
xlim([0.5 7.5])
xlabel('Day')
ylabel('Place cell ratio')


Cell_Cetegory_temp = ( cellfun(@(x) table2array(x), Cell_Cetegory_Day, 'UniformOutput', false) );
A = cell2mat(cellfun(@(x) sum(x(1:N_NegAnimal,4),'omitmissing'), Cell_Cetegory_temp, 'UniformOutput', false));
B = cell2mat(cellfun(@(x) sum(x(1:N_NegAnimal,3),'omitmissing'), Cell_Cetegory_temp, 'UniformOutput', false));
Neg = A./(A+B) * 100;
A = cell2mat(cellfun(@(x) sum(x(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,4),'omitmissing'), Cell_Cetegory_temp, 'UniformOutput', false));
B = cell2mat(cellfun(@(x) sum(x(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,3),'omitmissing'), Cell_Cetegory_temp, 'UniformOutput', false));
Pos = A./(A+B) * 100;
p = plot(Neg, '-bo', 'LineWidth', 3);
p = plot(Pos, '-ro', 'LineWidth', 3);


s = strcat("CellCategory_AcrossDay_AllCell.jpg");
exportgraphics(fig, fullfile(OutDir, s));



%% CellCategory_AcrossDay (per animal)

Cell_Cetegory_temp = cellfun(@(x) table2array(x), Cell_Cetegory_Day, 'UniformOutput', false);

nDay = numel(Cell_Cetegory_temp);

Neg_ratio = NaN(N_NegAnimal, nDay);
Pos_ratio = NaN(N_PosAnimal, nDay);

for d = 1:nDay

    data = Cell_Cetegory_temp{d};

    % ----- Negative animals -----
    place = data(1:N_NegAnimal,4);
    nonplace = data(1:N_NegAnimal,3);

    Neg_ratio(:,d) = place ./ (place + nonplace) * 100;


    % ----- Positive animals -----
    place = data(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,4);
    nonplace = data(N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,3);

    Pos_ratio(:,d) = place ./ (place + nonplace) * 100;

end


Neg_mean = mean(Neg_ratio,1,'omitnan');
Pos_mean = mean(Pos_ratio,1,'omitnan');

Neg_sem = std(Neg_ratio,0,1,'omitnan') ./ sqrt(sum(~isnan(Neg_ratio),1));
Pos_sem = std(Pos_ratio,0,1,'omitnan') ./ sqrt(sum(~isnan(Pos_ratio),1));

fig = figure;
set(fig, 'Position', [965   240   498   345]);
hold on

x = 1:nDay;

errorbar(x, Neg_mean, Neg_sem, '-bo', 'LineWidth',3)
errorbar(x, Pos_mean, Pos_sem, '-ro', 'LineWidth',3)

ylim([0 100])
xlim([0.5 7.5])
xlabel('Day')
ylabel('Place cell ratio')


x = 1:size(Neg_ratio,2);
% ----- Negative animals -----
for i = 1:N_NegAnimal
    y = Neg_ratio(i,:);
    valid = ~isnan(y);

    if any(valid)
        p = plot(x(valid), y(valid), '-o');
        p.Color = [0 0 1 0.2];
    end

end

x = 1:size(Pos_ratio,2);
% ----- Positive animals -----
for i = 1:N_PosAnimal
    y = Pos_ratio(i,:);
    valid = ~isnan(y);
    if any(valid)
        p = plot(x(valid), y(valid), '-o');
        p.Color = [1 0 0 0.2];
    end

end

s = strcat("CellCategory_AcrossDay_PerAnimal.jpg");
exportgraphics(fig, fullfile(OutDir, s));


%% FiringRate
close all

fieldName = 'FR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = ...
    func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);


for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

figure(fig2)
ylim([0 0.04])
yticks([0 0.02 0.04])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([0 0.1])
xticks([0 0.05 0.1])
yticklabels([])
xticklabels([])

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'SpInfo_meanRate';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

%
figure(fig2)
ylim([0 0.04])
yticks([0 0.02 0.04])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])

figure(fig3)
xlim([0 0.2])
xticks([0 0.1 0.2])
yticklabels([])
xticklabels([])

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%% spatial info
close all
fieldName = 'SpInfo_Zscore';
fieldNameShuf = 'SpInfo_shuf';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldNameShuf, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

%
figure(fig2)
ylim([0 6])
yticks([0 3 6])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([-2 14])
xticks([-2 0 7 14])
yticklabels([])
xticklabels([])

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));

s = strcat(fieldName, " cdf.pdf");
exportgraphics(fig3, fullfile(OutDir_panel, s));

%%

close all
fieldName = 'SpInfo_persec_Zscore';
fieldNameShuf = 'SpInfo_sec_shuf';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldNameShuf, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end


%
figure(fig2)
ylim([0 6])
yticks([0 3 6])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([-2 14])
xticks([-2 0 7 14])
yticklabels([])
xticklabels([])

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));



%%

close all
fieldName = 'SpInfo_spatialSelectivity';
fieldNameShuf = 'SpInfo_spatialSelectivity_shuf';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf_2(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldNameShuf, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

%
figure(fig2)
ylim([0 .8])
yticks([0 .4 .8])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([0 1])
xticks([0 .5 1])
yticklabels([])
xticklabels([])

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.pdf");
exportgraphics(fig3, fullfile(OutDir_panel, s));




%% Area
close all
fieldName = 'Area';
fieldName2 = 'BorderScore';
fieldThr = 0.5;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig2_2, fig3, fig3_2] = ...
    func_ShowCDFandBox_withF2(CellListT_all_animals, CellList_category_all, FieldStats_table_animal, fieldName, fieldName2, fieldThr, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
j = 1;
for i = 1:11
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end
j = 2;
for i = 1:11
    s = strcat(fieldName, "_TwoBars ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end

%
figure(fig2)
ylim([0 2000])
yticks([0 1000 2000])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])

figure(fig2_2)
ylim([0 2000])
yticks([0 1000 2000])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 5.3])

figure(fig3)
xlim([0 7200])
xticks([0 3600 7200])
yticklabels([])
xticklabels([])

figure(fig3_2)
xlim([0 7200])
xticks([0 3600 7200])
yticklabels([])
xticklabels([])


s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, " cdf.pdf");
exportgraphics(fig3, fullfile(OutDir_panel, 'Field', s));

s = strcat(fieldName, "_border bar.emf");
exportgraphics(fig2_2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, "_border cdf.pdf");
exportgraphics(fig3_2, fullfile(OutDir_panel, 'Field', s));



%% other figures
%%
close all
% fieldName = 'Coherence_Z';
fieldName = 'Coherence_R';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);

for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

%
figure(fig2)
ylim([0 0.6])
yticks([0 0.3 0.6])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([0 0.8])
xticks([0 0.4 0.8])
yticklabels([])
xticklabels([])

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));



%%
close all
fieldName = 'Field_InfieldFR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

%
figure(fig2)
ylim([0 0.08])
yticks([0 0.04 0.08])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([0 0.24])
xticks([0 0.12 0.24])
yticklabels([])
xticklabels([])


s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));



%%
fieldName = 'Field_OutfieldFR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

%
figure(fig2)
ylim([0 0.02])
yticks([0 0.01 0.02])

yticklabels([])
xticklabels([])
title('')
xlim([0.2 2.8])


figure(fig3)
xlim([0 0.08])
xticks([0 0.04 0.08])
yticklabels([])
xticklabels([])

%
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));





%% field size
close all
fieldName = 'Field_SizeAll';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);

for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Field_SizeAll_ratioArena';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Field_SizeAll_ratioActive';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Field_N';
showhist = 1;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, showhist);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));



%%
close all

fieldName = 'Field_MaxBorderScore';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Field_PeakFR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Field_InfieldFR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


fieldName = 'Field_OutfieldFR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Field_ACG_distCPNN';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Field_ACG_angleDeg';
showhist = 1;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, showhist);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%% stability

close all
fieldName = 'Stab_Pearson_Z';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);

for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Stab_XCorr_Shift';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);

for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Rel_boot_R_med';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));



%% plot with shuffle data


%
close all
fieldName = 'Stab_Pearson_R';
fieldNameShuf = 'Stab_Pearson_R_shuf';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldNameShuf, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Rel_boot_R_med';
fieldNameShuf = 'Rel_BootR_shuf';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldNameShuf, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%% plot with border score categorization

%
close all

fieldName = 'MinD';
fieldName2 = 'BorderScore';
fieldThr = 0.5;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig2_2, fig3, fig3_2] =  ...
    func_ShowCDFandBox_withF2(CellListT_all_animals, CellList_category_all, FieldStats_table_animal, fieldName, fieldName2, fieldThr, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
j = 1;
for i = 1:11
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end
j = 2;
for i = 1:11
    s = strcat(fieldName, "_TwoBars ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, 'Field', s));

s = strcat(fieldName, "_border bar.emf");
exportgraphics(fig2_2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, "_border cdf.emf");
exportgraphics(fig3_2, fullfile(OutDir_panel, 'Field', s));

%
close all

fieldName = 'PeakFR';
fieldName2 = 'BorderScore';
fieldThr = 0.5;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig2_2, fig3, fig3_2] = ...
    func_ShowCDFandBox_withF2(CellListT_all_animals, CellList_category_all, FieldStats_table_animal, fieldName, fieldName2, fieldThr, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);


j = 1;
for i = 1:11
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end
j = 2;
for i = 1:11
    s = strcat(fieldName, "_TwoBars ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, 'Field', s));

s = strcat(fieldName, "_border bar.emf");
exportgraphics(fig2_2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, "_border cdf.emf");
exportgraphics(fig3_2, fullfile(OutDir_panel, 'Field', s));


%%
close all

fieldName = 'BorderScore';
fieldName2 = 'BorderScore';
fieldThr = 0.5;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig2_2, fig3, fig3_2] = ...
    func_ShowCDFandBox_withF2(CellListT_all_animals, CellList_category_all, FieldStats_table_animal, fieldName, fieldName2, fieldThr, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);

j = 1;
for i = 1:11
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end
j = 2;
for i = 1:11
    s = strcat(fieldName, "_TwoBars ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, 'Field', s));

s = strcat(fieldName, "_border bar.emf");
exportgraphics(fig2_2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, "_border cdf.emf");
exportgraphics(fig3_2, fullfile(OutDir_panel, 'Field', s));


%%
close all

fieldName = 'NumPeak';
fieldName2 = 'BorderScore';
fieldThr = 0.5;
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig2_2, fig3, fig3_2] = ...
    func_ShowCDFandBox_withF2(CellListT_all_animals, CellList_category_all, FieldStats_table_animal, fieldName, fieldName2, fieldThr, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);

j = 1;
for i = 1:11
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end
j = 2;
for i = 1:11
    s = strcat(fieldName, "_TwoBars ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i, j}, fullfile(OutDir, 'Field', s));
end

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, 'Field', s));

s = strcat(fieldName, "_border bar.emf");
exportgraphics(fig2_2, fullfile(OutDir_panel, 'Field', s));
s = strcat(fieldName, "_border cdf.emf");
exportgraphics(fig3_2, fullfile(OutDir_panel, 'Field', s));


%% speed modulation
[status,msg,msgID] = mkdir(fullfile(OutDir, 'Speed'));

close all
fieldName = 'Speed_SpearmanR';
fieldName2 = 'Speed_Spearman_shuffP';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_Pvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));

%
close all
fieldName = 'Speed_SpearmanR';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_R';
fieldName2 = 'Speed_R_shuff_p';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_Pvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));

%
close all
fieldName = 'Speed_R';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_Info_spike';
fieldName2 = 'Speed_Info_p';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_Pvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_Info_spike';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_Info_spike';
fieldName2 = 'Speed_Spearman_shuffP';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_SpearmanPvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_Info_spike_z';
fieldName2 = 'Speed_Info_p';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_Pvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_Info_spike_z';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldName, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


close all
fieldName = 'Speed_Info_spike_z';
fieldName2 = 'Speed_Spearman_shuffP';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_SpearmanPvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%
close all
fieldName = 'Speed_Stab_r';
fieldName2 = 'Speed_Spearman_shuffP';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, 0, fieldName2, 0.05);
for i = 1:length(fig)
    s = strcat(fieldName, "_Pvalue ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end
s = strcat(fieldName, " p bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " p cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));


%%
close all
fieldName = 'Rel_speed_boot_R_med';
fieldNameShuf = 'Rel_speed_BootR_shuf';
[fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
for i = 1:length(fig)
    s = strcat(fieldNameShuf, " ", sgst_fig{i}, ".jpg");
    exportgraphics(fig{i}, fullfile(OutDir, 'Speed', s));
end

s = strcat(fieldName, " bar.emf");
exportgraphics(fig2, fullfile(OutDir_panel, s));
s = strcat(fieldName, " cdf.emf");
exportgraphics(fig3, fullfile(OutDir_panel, s));



%% percent encoding speed
fieldThr = 0.05;

close all

fieldName = 'Speed_Info_p';
Ymax = 100;
[fig, sgst_fig, pval, fig2] = func_CellPercentage(CellListT_all_animals, CellList_category_all, fieldName, fieldThr, Ymax, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
s = strcat(sgst_fig, ".jpg");
exportgraphics(fig, fullfile(OutDir, 'Speed', s));
s = strcat(sgst_fig, "_IndAnimal.jpg");
exportgraphics(fig2, fullfile(OutDir, 'Speed', s));


fieldName = 'Speed_Spearman_shuffP';
Ymax = 40;
[fig, sgst_fig, pval, fig2] = func_CellPercentage(CellListT_all_animals, CellList_category_all, fieldName, fieldThr, Ymax, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
s = strcat(sgst_fig, ".jpg");
exportgraphics(fig, fullfile(OutDir, 'Speed', s));
s = strcat(sgst_fig, "_IndAnimal.jpg");
exportgraphics(fig2, fullfile(OutDir, 'Speed', s));


fieldName = 'Speed_R_shuff_p';
Ymax = 40;
[fig, sgst_fig, pval, fig2] = func_CellPercentage(CellListT_all_animals, CellList_category_all, fieldName, fieldThr, Ymax, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN);
s = strcat(sgst_fig, ".jpg");
exportgraphics(fig, fullfile(OutDir, 'Speed', s));
s = strcat(sgst_fig, "_IndAnimal.jpg");
exportgraphics(fig2, fullfile(OutDir, 'Speed', s));

%%

close all


%% save variables for Naoki

s_save = strcat('PlaceCellInfo_260218.mat');
save(s_save,...
    'mice_str', 'CellListT_all_animals',...
    'CellList_category_all', 'FieldStats_table_animal',...
    '-v7.3');


%% functions
%% Generate Cell List
function [CellFieldList_table, CellFieldList_shuf_table, CellSpeedModList_table] = func_GenerateCellList(LoadData, Boot, mice_ind, Day)
%% Load Basics
SpikeN_All  = cellfun(@length, LoadData{mice_ind, Day}.BasicAnalysis.SpkDetectionResults_save.Peaks');
CellInd     = (1:length(SpikeN_All))';
Animal      = ones(length(SpikeN_All),1) * mice_ind;

runthr = 2;
SpikeN_Run_table    = LoadData{mice_ind, Day}.BasicAnalysis.SpkDetectionResults_save.SpkTimeAligned_table;
SpikeN_Run_N        = cellfun(@(x) sum(x.V_1secMean > runthr), SpikeN_Run_table)';

FR = SpikeN_All / LoadData{mice_ind, Day}.BasicAnalysis.TotalRecordingTime;

Rel_boot_p      = Boot{mice_ind, Day}.p_rank_ratemap_cell;
Rel_boot_delta  = Boot{mice_ind, Day}.Cliff_delta_ratemap_cell;
Rel_boot_R_med  = cellfun(@(x) median(x, 'omitmissing'), Boot{mice_ind, Day}.HalvesCorr_Pearson_Boot_cell);
Rel_boot_R_CI95 = Boot{mice_ind, Day}.CI_rateMap_cell(:,1);

SpInfo      = cellfun(@(x) x.info_obs, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);
SpInfo_p    = cellfun(@(x) x.pval_spinfo, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);

SpInfo_Zscore      = cellfun(@(x) x.SpInfo_z, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);
SpInfo_Normal      = cellfun(@(x) x.SpInfo_n, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);

SpInfo_meanRate             = cellfun(@(x) x.meanRate, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_cell);
SpInfo_persec               = cellfun(@(x) x.info_persec, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_cell);
SpInfo_persec_p             = cellfun(@(x) x.pval_spsecinfo, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);
SpInfo_spatialSparsity      = cellfun(@(x) x.spatialSparsity, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_cell);
SpInfo_spatialSelectivity      = cellfun(@(x) x.spatialSelectivity, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_cell);
SpInfo_spatialSparsity_p    = cellfun(@(x) x.pval_sparsity, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);
SpInfo_persec_Zscore_temp               = cellfun(@(x) x.info_persec, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_cell);
SpInfo_persec_Zscore_temp_mean          = cellfun(@(x) mean(x.info_sec_shuf, 'omitmissing'), LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);
SpInfo_persec_Zscore_temp_std          = cellfun(@(x) std(x.info_sec_shuf, 'omitmissing'), LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell);
SpInfo_persec_Zscore = (SpInfo_persec_Zscore_temp - SpInfo_persec_Zscore_temp_mean) ./ SpInfo_persec_Zscore_temp_std;

Coherence_R = cellfun(@(x) x.SpCorr, LoadData{mice_ind, Day}.BasicAnalysis.Coherence_cell);
Coherence_Z = cellfun(@(x) x.Z, LoadData{mice_ind, Day}.BasicAnalysis.Coherence_cell);

Field_N                = cellfun(@(x) x.FieldNum, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);
Field_MaxBorderScore   = cellfun(@(x) x.MaxBorderScore, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);
Field_PeakFR           = cellfun(@(x) x.PeakFR, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);
Field_InfieldFR        = cellfun(@(x) x.infield_meanfiring, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);
Field_OutfieldFR       = cellfun(@(x) x.outfield_meanfiring, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);

ArenaArea = LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell{1}.ArenaArea;
ActiveArea = LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell{1}.ActiveArea;
Field_SizeAll              = cellfun(@(x) x.FieldSize_All, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);
Field_SizeAll_ratioArena   = Field_SizeAll/ArenaArea *100;
Field_SizeAll_ratioActive  = Field_SizeAll/ActiveArea *100;
Field_SizePeak             = cellfun(@(x) x.FieldSize_Peak, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell);
Field_SizePeak_ratioArena  = Field_SizePeak/ArenaArea *100;
Field_SizePeak_ratioActive = Field_SizePeak/ActiveArea *100;

Field_ACG_distCPNN  = cellfun(@(x) x.distCPNN, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_ACGanalysis_cell);
Field_ACG_angleDeg  = cellfun(@(x) x.angleDeg, LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_ACGanalysis_cell);
Field_ACG_angleDeg(Field_ACG_angleDeg>179.9) = NaN;

speedThrID = 3; %2cm/s
Stability_FR = cell2mat((cellfun(@(x) x.FiringRate_halve, LoadData{mice_ind, Day}.BasicAnalysis.HalveCompareField_stability_cell, 'UniformOutput', false))')';
Stab_FR_early = Stability_FR(:,1);
Stab_FR_later = Stability_FR(:,2);
Stab_Pearson_R = cell2mat(cellfun(@(x) x.HalvesCorr_Pearson_sp(speedThrID,1),   LoadData{mice_ind, Day}.BasicAnalysis.HalveCompareField_stability_cell, 'UniformOutput', false));
Stab_Pearson_Z = cell2mat(cellfun(@(x) x.HalvesCorr_Pearson_z_sp(speedThrID,1), LoadData{mice_ind, Day}.BasicAnalysis.HalveCompareField_stability_cell, 'UniformOutput', false));
Stab_Pearson_p = cell2mat(cellfun(@(x) x.HalvesCorr_Pearson_p_sp(speedThrID,1), LoadData{mice_ind, Day}.BasicAnalysis.HalveCompareField_stability_cell, 'UniformOutput', false));
Stab_XCorr_Max = cell2mat(cellfun(@(x) x.XCor_max_sp(speedThrID,1),   LoadData{mice_ind, Day}.BasicAnalysis.HalveCompareField_stability_cell, 'UniformOutput', false));
Stab_XCorr_Shift = cell2mat(cellfun(@(x) x.XCor_shiftDist_sp(speedThrID,1), LoadData{mice_ind, Day}.BasicAnalysis.HalveCompareField_stability_cell, 'UniformOutput', false));

Stab_shuff = cellfun(@(x) x.Field.shuf,   LoadData{mice_ind, Day}.BasicAnalysis.Halves_shuf_stab_cell);
Stab_Pearson_R_shufP = vertcat(Stab_shuff.HalvesCorr_Pearson_shuf_pval);


Rel_speed_boot_R_med      = cellfun(@(x) median(x, 'omitmissing'), Boot{mice_ind, Day}.HalvesCorr_speed_Pearson_Boot_cell);
Rel_speed_boot_p          = Boot{mice_ind, Day}.p_rank_speed_cell;


% summarize
CellList_array = [Animal, CellInd, SpikeN_All, SpikeN_Run_N, Rel_boot_p, Rel_boot_delta, Rel_boot_R_med, Rel_boot_R_CI95, FR, SpInfo, SpInfo_p,...
    SpInfo_Zscore, SpInfo_Normal, SpInfo_meanRate, SpInfo_persec, SpInfo_persec_p, SpInfo_spatialSparsity, SpInfo_spatialSelectivity, SpInfo_spatialSparsity_p,...
    SpInfo_persec_Zscore, ...
    Coherence_R, Coherence_Z, ...
    Field_N, Field_MaxBorderScore, Field_PeakFR, Field_InfieldFR, Field_OutfieldFR, ...
    Field_SizeAll, Field_SizeAll_ratioArena, Field_SizeAll_ratioActive, Field_SizePeak, Field_SizePeak_ratioArena, Field_SizePeak_ratioActive, ...
    Field_ACG_distCPNN, Field_ACG_angleDeg,...
    Stab_FR_early, Stab_FR_later, Stab_Pearson_R, Stab_Pearson_Z, Stab_Pearson_p, Stab_XCorr_Max, Stab_XCorr_Shift, Stab_Pearson_R_shufP,...
    Rel_speed_boot_R_med, Rel_speed_boot_p];
varNames = {'Animal', 'CellInd', 'SpikeN_All', 'SpikeN_Run_N', 'Rel_boot_p', 'Rel_boot_delta', 'Rel_boot_R_med', 'Rel_boot_R_CI95', 'FR', 'SpInfo', 'SpInfo_p',...
    'SpInfo_Zscore', 'SpInfo_Normal', 'SpInfo_meanRate', 'SpInfo_persec', 'SpInfo_persec_p', 'SpInfo_spatialSparsity', 'SpInfo_spatialSelectivity', 'SpInfo_spatialSparsity_p',...
    'SpInfo_persec_Zscore', ...
    'Coherence_R', 'Coherence_Z', ...
    'Field_N', 'Field_MaxBorderScore', 'Field_PeakFR', 'Field_InfieldFR', 'Field_OutfieldFR', ...
    'Field_SizeAll', 'Field_SizeAll_ratioArena', 'Field_SizeAll_ratioActive', 'Field_SizePeak', 'Field_SizePeak_ratioArena', 'Field_SizePeak_ratioActive', ...
    'Field_ACG_distCPNN', 'Field_ACG_angleDeg',...
    'Stab_FR_early', 'Stab_FR_later', 'Stab_Pearson_R', 'Stab_Pearson_Z', 'Stab_Pearson_p', 'Stab_XCorr_Max', 'Stab_XCorr_Shift', 'Stab_Pearson_R_shufP',...
    'Rel_speed_boot_R_med', 'Rel_speed_boot_p'};

CellFieldList_table = array2table(CellList_array, 'VariableNames', varNames);

%% Load shuf data

SpInfo_shuf      = cellfun(@(x) x.info_shuf, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell, 'UniformOutput',false);
SpInfo_sec_shuf      = cellfun(@(x) x.info_sec_shuf, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell, 'UniformOutput',false);
SpInfo_spatialSparsity_shuf      = cellfun(@(x) x.spatialSparsity_shuf, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell, 'UniformOutput',false);
SpInfo_spatialSelectivity_shuf      = cellfun(@(x) 1 - x.spatialSparsity_shuf, LoadData{mice_ind, Day}.BasicAnalysis.SpInfo_shuf_cell, 'UniformOutput',false);


Stab_Pearson_R_shuf = cellfun(@(x) x.Field.shuf.HalvesCorr_Pearson_shuf,   LoadData{mice_ind, Day}.BasicAnalysis.Halves_shuf_stab_cell, 'UniformOutput', false);


Rel_BootR           = LoadData{mice_ind, Day}.BOOT.HalvesCorr_Pearson_Boot_cell;
Rel_BootR_shuf      = LoadData{mice_ind, Day}.BOOT.HalvesCorr_Pearson_Boot_cell_shuf;


Rel_speed_BootR           = LoadData{mice_ind, Day}.BOOT.HalvesCorr_speed_Pearson_Boot_cell;
Rel_speed_BootR_shuf      = LoadData{mice_ind, Day}.BOOT.HalvesCorr_speed_Pearson_Boot_cell_speed;


CellList_array = [num2cell(Animal), num2cell(CellInd), SpInfo_shuf, SpInfo_sec_shuf, SpInfo_spatialSelectivity_shuf, Stab_Pearson_R_shuf,...
    Rel_BootR, Rel_BootR_shuf, ...
    Rel_speed_BootR, Rel_speed_BootR_shuf];
varNames = {'Animal', 'CellInd', 'SpInfo_shuf', 'SpInfo_sec_shuf', 'SpInfo_spatialSelectivity_shuf', 'Stab_Pearson_R_shuf',...
    'Rel_BootR', 'Rel_BootR_shuf', ...
    'Rel_speed_BootR', 'Rel_speed_BootR_shuf'};

CellFieldList_shuf_table = array2table(CellList_array, 'VariableNames', varNames);

%% Load SpeedMod
Speed_R         = cellfun(@(x) x.corr_r, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_R_p       = cellfun(@(x) x.corr_p, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_R_shuff_p         = cellfun(@(x) x.shuffle_p, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);

Speed_SpearmanR         = cellfun(@(x) x.corr_SpearmanR, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_SpearmanP         = cellfun(@(x) x.corr_SpearmanP, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_Spearman_shuffP   = cellfun(@(x) x.p_shuffle_Spearman, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);

Speed_Info_spike        = cellfun(@(x) x.speedInfo.I_spike, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_Info_sec          = cellfun(@(x) x.speedInfo.I_sec, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_Info_p            = cellfun(@(x) x.speedInfo.I_spikeShuf_p, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);

meanTemp = cellfun(@(x) mean(x.speedInfo.I_spikeShuf, 'omitmissing'), LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
stdTemp  = cellfun(@(x) std(x.speedInfo.I_spikeShuf, 'omitmissing'), LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_Info_spike_z          = (Speed_Info_spike - meanTemp)./stdTemp;
Speed_Info_spike_normal     = Speed_Info_spike./meanTemp;
meanTemp = cellfun(@(x) mean(x.speedInfo.I_secShuf, 'omitmissing'), LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
stdTemp  = cellfun(@(x) std(x.speedInfo.I_secShuf, 'omitmissing'), LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_Info_sec_z          = (Speed_Info_sec - meanTemp)./stdTemp;
Speed_Info_sec_normal     = Speed_Info_sec./meanTemp;

Speed_slope = cellfun(@(x) x.regression_basics.slope, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);
Speed_LogisticSlope = cellfun(@(x) x.LogisticSlopes, LoadData{mice_ind, Day}.BasicAnalysis.speedModResult_cell);

Speed_Stab_r = cellfun(@(x) x.Speed.real.r_Pearson, LoadData{mice_ind, Day}.BasicAnalysis.Halves_shuf_stab_cell);
Speed_Stab_r_p = cellfun(@(x) x.Speed.real.p_Pearson, LoadData{mice_ind, Day}.BasicAnalysis.Halves_shuf_stab_cell);
Speed_Stab_logiSlopDiff = cellfun(@(x) x.Speed.real.stability_logistic, LoadData{mice_ind, Day}.BasicAnalysis.Halves_shuf_stab_cell);
Speed_Stab_r_shuffP = cellfun(@(x) x.Speed.shuf.pval_Pearson, LoadData{mice_ind, Day}.BasicAnalysis.Halves_shuf_stab_cell);


% summarize
CellSpeedModList_array = [CellInd, Speed_R, Speed_R_p, Speed_R_shuff_p, Speed_SpearmanR, Speed_SpearmanP, Speed_Spearman_shuffP,...
    Speed_Info_spike, Speed_Info_sec, Speed_Info_p, Speed_Info_spike_z, Speed_Info_spike_normal, Speed_Info_sec_z, Speed_Info_sec_normal, ...
    Speed_slope, Speed_LogisticSlope, Speed_Stab_r, Speed_Stab_r_p, Speed_Stab_logiSlopDiff, Speed_Stab_r_shuffP, ...
    ];
varNames = {'CellInd', 'Speed_R', 'Speed_R_p', 'Speed_R_shuff_p', 'Speed_SpearmanR', 'Speed_SpearmanP', 'Speed_Spearman_shuffP',...
    'Speed_Info_spike', 'Speed_Info_sec', 'Speed_Info_p', 'Speed_Info_spike_z', 'Speed_Info_spike_normal', 'Speed_Info_sec_z', 'Speed_Info_sec_normal', ...
    'Speed_slope', 'Speed_LogisticSlope', 'Speed_Stab_r', 'Speed_Stab_r_p', 'Speed_Stab_logiSlopDiff', 'Speed_Stab_r_shuffP', ...
    };


CellSpeedModList_table = array2table(CellSpeedModList_array, 'VariableNames', varNames);

end

%% Generate Field List
function FieldStats_table = func_GenerateFieldList(LoadData, mice_ind, Day, CellList_table)
data = LoadData{mice_ind, Day}.FieldAnalysisRedo.FRD_FieldAnalysis_Redo_cell;

% total number of fields
FieldN_AllCell = CellList_table.Field_N;
FieldN_AllCell(isnan(FieldN_AllCell)) = 1;
Nf = sum(FieldN_AllCell);

% Preallocate
FieldStats_Ind              = nan(Nf,1);
FieldStats_Area             = nan(Nf,1);
FieldStats_MinD             = nan(Nf,1);
FieldStats_Circularity      = nan(Nf,1);
FieldStats_PeakFR           = nan(Nf,1);
FieldStats_NumPeak          = nan(Nf,1);
FieldStats_BorderScore      = nan(Nf,1);
FieldStats_BorderDist       = nan(Nf,1);

count = 0;

for i = 1:numel(data)
    cid = data{i}.CellID;
    s   = data{i}.stats;

    % is stats is empty or 0x1 struct → continue
    if isempty(s) || (iscell(s) && isscalar(s) && isempty(s{1}))
        count = count + 1;
        FieldStats_Ind(count)  = cid;
        continue
    end

    % unify struct arrays to cells
    if isstruct(s)
        statsArr = num2cell(s);
    else
        statsArr = s;
    end

    % extract each fields
    for j = 1:numel(statsArr)
        fs = statsArr{j};
        if isstruct(fs)
            for k = 1:numel(fs)
                count = count + 1;
                if isfield(fs(k),'Area') && ~isempty(fs(k).Area)
                    FieldStats_Area(count) =        fs(k).Area_real;

                    if numel(fs) > 1
                        Pos = fs.Centroid;
                        PosAll = vertcat(fs.Centroid);
                        D = pdist2(Pos, PosAll);
                        MinD = min(D(D~=0));
                        FieldStats_MinD(count) = MinD;
                    end
                    FieldStats_Circularity(count) = fs(k).Circularity;
                    FieldStats_PeakFR(count) = fs(k).MaxIntensity;
                    FieldStats_NumPeak(count) = size(fs(k).LocalPeakI,1);
                    FieldStats_BorderScore(count) = fs(k).BorderScore_PolygonArena;
                    FieldStats_BorderDist(count) = fs(k).BorderScore_PolygonArena_Dist;

                else
                end
                FieldStats_Ind(count) = cid;
            end
        else
            % just in case of error
            count = count + 1;
            FieldStats_Ind(count)  = cid;
        end
    end
end

% summarize
Animal = ones(Nf,1) * mice_ind;
FieldStats = [Animal, FieldStats_Ind, FieldStats_Area, FieldStats_MinD, FieldStats_Circularity, ...
    FieldStats_PeakFR, FieldStats_NumPeak, FieldStats_BorderScore, FieldStats_BorderDist];
varNames = {'Animal', 'CellInd', 'Area', 'MinD', 'Circularity', ...
    'PeakFR', 'NumPeak', 'BorderScore', 'BorderDist'};
FieldStats_table = array2table(FieldStats, 'VariableNames',varNames);

end

%% show cdf plot and box plot
%% func_ShowCDFandBox
function [fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] = func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, fieldName, N_NegAnimal, N_PosAnimal,ind_NegNonNaN,ind_PosNonNaN, varargin)
%   1. plot CDF (default)
%       func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, 'SpatialInfo')
%
%   2. show histogram
%       func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, 'SpatialInfo', 1)
%
%   3. use condition filtering with fieldName2
%       func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, 'SpatialInfo', 1, 'SpeedModulation', 0.3)


showhist = 0;
useField2 = false;

if nargin >= 8 && ~isempty(varargin{1})
    showhist = varargin{1};
end

if nargin >= 9 && ~isempty(varargin{2})
    fieldName2 = varargin{2};
    useField2 = true;
    fieldThr = varargin{3};
end


%%
fig = cell(11,1);
Neg_allDay = cell(7,5);
Pos_allDay = cell(7,5);

sgst_fig    = cell(11,1);
p_ks        = nan(10,5);
ks2stat     = nan(10,5);
p_rank      = nan(10,5);
stats_rank  = cell(10,5);
e_delta     = nan(10,5);

for Day = 1:11
    %%

    % fig{Day} = figure('Position', [-1584, 404, 1486, 440], 'Visible','off');
    fig{Day} = figure('Position', [-1584, 404, 1486, 440], 'Visible','on');
    figrow = 2;
    figcol = 6;

    % define fields
    fields = { ...
        'CellListT_all_animals', 'Cell_unstableRatemap', ...
        'Cell_reliableRatemap', 'Cell_reliableRatemap_NotPlaceCell', ...
        'Cell_reliableRatemap_PlaceCell'};

    % fig indices
    fig_indices = [1, 3, 5, 7, 9];

    % show plots
    for k = 1:numel(fields)
        if Day<=7
            fieldname = fields{k};
            if k == 1
                A = CellListT_all_animals(:, Day);
                AA = A;
            else
                A = CellList_category_all(:, Day);

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end
            end

            % grouping Negative and Positive
            Neg = vertcat(AA{1:N_NegAnimal,:});
            Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

            Neg_allDay{Day,k} = Neg;
            Pos_allDay{Day,k} = Pos;

        end

        if Day == 8
            Neg = vertcat(Neg_allDay{:, k});
            Pos = vertcat(Pos_allDay{:, k});
        elseif Day == 9
            Neg = vertcat(Neg_allDay{1:3, k});
            Pos = vertcat(Pos_allDay{1:3, k});
        elseif Day == 10
            Neg = vertcat(Neg_allDay{4:7, k});
            Pos = vertcat(Pos_allDay{4:7, k});
        end

        %
        if Day == 11 % three successful days

            if k == 1
                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];

                AA = CellListT_all_animals;
                AA(~ind_temp) = {[]};

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            else

                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];
                A = CellList_category_all;
                A(~ind_temp) = {[]};


                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

                Neg_3trials = Neg;
                Pos_3trials = Pos;
                sgst = strcat("three successful trials");

            end

        end
        %%

        % === field data ===
        if isempty(Neg)
            Data_A = [];
        else
            Data_A = Neg.(fieldName);
        end

        if isempty(Pos)
            Data_B = [];
        else
            Data_B = Pos.(fieldName);
        end

        % === condition filtering ===
        if useField2
            if isempty(Neg)
                Data_A = [];
            else
                condA = Neg.(fieldName2) < fieldThr;
                Data_A = Data_A(condA);
            end

            if isempty(Pos)
                Data_B = [];
            else
                condB = Pos.(fieldName2) < fieldThr;
                Data_B = Data_B(condB);
            end
        end


        % === show plots ===
        fig_i = fig_indices(k);
        [p_ks(Day, k), ks2stat(Day, k), p_rank(Day, k), stats_rank{Day, k}, e_delta(Day, k)] = ...
            func_plotCdgAndBox(Data_A, Data_B, figrow, figcol, fig_i, fieldName, showhist);

        if Day<=7
            sgst = strcat("Day", num2str(Day));
        elseif Day == 8
            sgst = strcat("Day All");
        elseif Day == 9
            sgst = strcat("Day 1-3");
        elseif Day == 10
            sgst = strcat("Day 4-7");
        end
        sgtitle(sgst)
        sgst_fig{Day} = sgst;
    end


    titles = {'All cell', 'Unstable cell', 'Reliable cell', 'Non-place cell', 'Place cell'};
    for i = 1:numel(titles)
        for j = 1:2
            idx = (i-1)*2 + j;
            subplot(2,6,idx);
            ax = gca;
            ax.Title.String = {titles{i}, ax.Title.String};
        end
    end

    for i = 1:5
        subplot(2,6, i*2)
        ylabel(fieldName, 'Interpreter', 'none')
    end


end



%% For fig panel
Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

% Day1-3
Neg = Neg_3trials;
Pos = Pos_3trials;

% === field data ===
if isempty(Neg)
    Data_A = [];
else
    Data_A = Neg.(fieldName);
end
if isempty(Pos)
    Data_B = [];
else
    Data_B = Pos.(fieldName);
end

% === condition filtering ===
if useField2
    condA = Neg.(fieldName2) < fieldThr;
    Data_A = Data_A(condA);
    condB = Pos.(fieldName2) < fieldThr;
    Data_B = Data_B(condB);
end


% close all
fig2 = figure;
% get(gcf,'Position')
set(fig2, 'Position', [1120         708         120         170])

hold on


A = Data_A;
B = Data_B;

mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

hold on
bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
ax = gca;
ax.Layer = 'top';

errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);
errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);



xticks(1:2);
xticklabels({'Neg', 'Pos'})


[h,p,ci,stats] = ttest2(A, B);

title('')
xlim([0.2 2.8])


disp(mA)
disp(mB)
disp(semA)
disp(semB)


sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Control: %.3f ± %.3f\n', mA, sdA);
fprintf('Casp: %.3f ± %.3f\n', mB, sdB);

%%

fig3 = figure;
% get(gcf,'Position')
set(fig3, 'Position', [1271         708         203         170])

hold on


hold on;
if isempty(Data_A)
    h1 = cdfplot(0);
else
    h1 = cdfplot(Data_A);
end
if isempty(Data_B)
    h2 = cdfplot(0);
else
    h2 = cdfplot(Data_B);
end
h1.Color = 'k';
h1.LineWidth = 2;
h2.Color = Ccasp;
h2.LineWidth = 2;

grid off;

try
    [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
    title(sprintf('K–S test p = %.4f', p_ks));
catch
    p_ks = 1; ks2stat = 0;
end

title('')
ylim([0 1])
yticks([0 0.5 1])
xlabel('')
ylabel('')



end


function [p_ks, ks2stat, p_rank, stats_rank, e_delta] = func_plotCdgAndBox(Data_A, Data_B, figrow, figcol, fig_i, fieldName, showhist)
%% cdf plot

if isempty(Data_A)
    Data_A = 0;
end
if isempty(Data_B)
    Data_B = 0;
end


if showhist == 0
    subplot(figrow, figcol, fig_i);
    hold on;
    if isempty(Data_A)
        h1 = cdfplot(0);
    else
        h1 = cdfplot(Data_A);
    end
    if isempty(Data_B)
        h2 = cdfplot(0);
    else
        h2 = cdfplot(Data_B);
    end
    h1.Color = 'b';
    h1.LineWidth = 1;
    h2.Color = 'r';
    h2.LineWidth = 1;
    lgd = legend({'Negative', 'Positive'},'Location','southeast');
    lgd.ItemTokenSize = [10, 8];

    xlabel(fieldName, 'Interpreter','none')
    ylabel('Cumulative probability');
    grid off;

    % --- Kolmogorov–Smirnov test ---
    if isempty(Data_A) == 0 && isempty(Data_B) == 0
        try
            [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
            title(sprintf('K–S test p = %.4f', p_ks));
        catch
            p_ks = 1; ks2stat = 0;
        end
    else
        p_ks = 1; ks2stat = 0;
    end
    title(sprintf('K–S test p = %.4f', p_ks));

    func_showstar_forCDF(p_ks);

elseif showhist == 1
    %%
    subplot(figrow, figcol, fig_i);
    hold on;

    Max = max([Data_A(:); Data_B(:)]);
    edges = -0.5:1:Max-0.5;
    if Max == 0
        edges = -0.5:0.5;
    end
    if length(edges)>20
        edges = 0:10:180;
    end
    binCenters = edges(1:end-1) + diff(edges)/2;
    [counts1,~] = histcounts(Data_A, edges);
    [counts2,~] = histcounts(Data_B, edges);

    counts1 = counts1/length(Data_A(~isnan(Data_A)));
    counts2 = counts2/length(Data_B(~isnan(Data_B)));
    dataMat = [counts1(:), counts2(:)];

    % figure;
    b = bar(binCenters, dataMat, 'grouped');
    b(1).FaceColor = [0 0 1];
    b(2).FaceColor = [1 0 0];
    b(1).FaceAlpha = 0.5;
    b(2).FaceAlpha = 0.5;
    xlabel(fieldName, 'Interpreter','none');
    ylabel('Probability');
    title('Side-by-side histogram (normalized to probability)');

    try
        xlim([binCenters(1), binCenters(end)] + 0.5)
        ylim([0, max([counts1(:); counts2(:)]) * 1.2])
    catch
    end

    box off;

    % --- Kolmogorov–Smirnov test ---
    if isempty(Data_A) == 0 && isempty(Data_B) == 0
        try
            [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
            title(sprintf('K–S test p = %.4f', p_ks));
        catch
            p_ks = 1; ks2stat = 0;
        end
    else
        p_ks = 1; ks2stat = 0;
    end

    title(sprintf('K–S test p = %.4f', p_ks));
    func_showstar_forCDF(p_ks);

end

%% box plot

subplot(figrow, figcol, fig_i+1)
hold on
b = boxplot([Data_A; Data_B], ...
    [repmat({'Negative'}, numel(Data_A), 1); ...
    repmat({'Positive'}, numel(Data_B), 1)], ...
    'Colors', [0 0 1; 1 0 0], 'Symbol', 'o', 'Whisker', 1.5);

ylabel(fieldName, 'Interpreter','none')
title('Spatial Information Comparison');
set(gca, 'Box', 'off', 'FontSize', 8);
xlim([0.2 2.8])

ax = gca;
Ymax = prctile([Data_A(:); Data_B(:)], 98);
ax.YLim = [ax.YLim(1), Ymax*1.4];



try
    [p_rank, ~, stats_rank] = ranksum(Data_A, Data_B);
    Effect = meanEffectSize(Data_A, Data_B, "Effect", "cliff" );
    e_delta = Effect.Effect;

    n1 = numel(Data_A);
    n2 = numel(Data_B);
    [p, ~, stats] = ranksum(Data_A, Data_B);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('All U1 = %.0f, U2 = %.0f, p = %.6g\n', U1, U2, p);
    stats_rank.U1 = U1;
    stats_rank.U2 = U2;
catch
    p_rank = 1; stats_rank = [];
    e_delta = 0;
end

func_showstar_box(p_rank, e_delta)

end

%
function func_showstar_box(p_val, e_delta)
ax = gca;
yLim = ax.YLim;
yStar = yLim(2) + range(yLim)*0.1;
ylim([yLim(1), yStar + range(yLim)*0.3]);

p = p_val;
if p < 0.05 && abs(e_delta) > 0.147
    plot([1 2], [yStar yStar], 'k-', 'LineWidth', 1);

    if abs(e_delta) > 0.474
        text(1.5, yStar + range(yLim)*0.05, '###', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    elseif abs(e_delta) <= 0.474 && abs(e_delta) > 0.33
        text(1.5, yStar + range(yLim)*0.05, '##', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    elseif abs(e_delta) <= 0.33 && abs(e_delta) > 0.147
        text(1.5, yStar + range(yLim)*0.05, '#', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    end

    % p_str = sprintf('p=%.3f', p_val);
    % e_str = sprintf('δ=%.3f', e_delta);
    % text(1.5, yStar + range(yLim)*0.20, {'ignorable'}, ...
    %     'FontSize', 10, 'HorizontalAlignment', 'center');

elseif p < 0.05 && abs(e_delta) <= 0.147
    text(1.5, yStar + range(yLim)*0.20, {'ignorable'}, ...
        'FontSize', 10, 'HorizontalAlignment', 'center');

elseif p>=0.05
    % p_str = sprintf('p=%.3f', p_val);
    % e_str = sprintf('δ=%.3f', e_delta);
    text(1.5, yStar + range(yLim)*0.20, {'n.s.'}, ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
end

end

%
function func_showstar_forCDF(p_val)
ax = gca;
yLim = ax.YLim;
yStar = yLim(2)*0.9;

p = p_val;
xLim = ax.XLim;
xpos = xLim(1) + range(xLim)*0.1;

% range(yLim);
if p < 0.05
    if p <0.01 && p>=0.001
        text(xpos, yStar + range(yLim)*0.05, '**', ...
            'FontSize', 8, 'HorizontalAlignment', 'center');
    elseif p<0.001
        text(xpos, yStar + range(yLim)*0.05, '***', ...
            'FontSize', 8, 'HorizontalAlignment', 'center');
    else
        text(xpos, yStar + range(yLim)*0.05, '*', ...
            'FontSize', 8, 'HorizontalAlignment', 'center');
    end

else

    text(xpos, yStar + range(yLim)*0.05, {'n.s.'}, ...
        'FontSize', 8, 'HorizontalAlignment', 'center');
end

end


%%
%% func_ShowCDFandBox with shuffle data
function [fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, varargin)
%   1. plot CDF (default)
%       func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, 'SpatialInfo')
%   2. show histogram
%       func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, 'SpatialInfo', 1)
%   3. use condition filtering with fieldName2
%       func_ShowCDFandBox(CellListT_all_animals, CellList_category_all, 'SpatialInfo', 1, 'SpeedModulation', 0.3)


showhist = 0;
useField2 = false;

if nargin >= 11 && ~isempty(varargin{1})
    showhist = varargin{1};
end


%%
fig = cell(11,1);
Neg_allDay = cell(7,5);
Pos_allDay = cell(7,5);
Neg_allDay_shuf = cell(7,5);
Pos_allDay_shuf = cell(7,5);

sgst_fig    = cell(11,1);
p_ks        = nan(10,5);
ks2stat     = nan(10,5);
p_rank      = nan(10,5);
stats_rank  = cell(10,5);
e_delta     = nan(10,5);

for Day = 1:11

    fig{Day} = figure('Position', [-1584, 404, 1486, 440]);
    figrow = 2;
    figcol = 6;

    % define fields
    fields = { ...
        'CellListT_all_animals', 'Cell_unstableRatemap', ...
        'Cell_reliableRatemap', 'Cell_reliableRatemap_NotPlaceCell', ...
        'Cell_reliableRatemap_PlaceCell'};

    % fig indices
    fig_indices = [1, 3, 5, 7, 9];

    % show plots
    for k = 1:numel(fields)
        Neg_shuf = [];
        Pos_shuf = [];


        if Day<=7
            fieldname = fields{k};
            if k == 1
                A = CellListT_all_animals(:, Day);
                AA = A;

                AA_shuf = CellListT_all_animals_shuf(:, Day);
            else
                A = CellList_category_all(:, Day);

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end

            end

            % grouping Negative and Positive
            Neg = vertcat(AA{1:N_NegAnimal,:});
            Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            Neg_allDay{Day,k} = Neg;
            Pos_allDay{Day,k} = Pos;

            Neg_shuf = vertcat(AA_shuf{1:N_NegAnimal,:});
            Pos_shuf = vertcat(AA_shuf{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            Neg_allDay_shuf{Day,k} = Neg_shuf;
            Pos_allDay_shuf{Day,k} = Pos_shuf;
        end

        if Day == 8
            Neg = vertcat(Neg_allDay{:, k});
            Pos = vertcat(Pos_allDay{:, k});
            % Neg_shuf = vertcat(Neg_allDay_shuf{:, k}); % very heavy to plot, use memory a lot
            % Pos_shuf = vertcat(Neg_allDay_shuf{:, k});
        elseif Day == 9
            Neg = vertcat(Neg_allDay{1:3, k});
            Pos = vertcat(Pos_allDay{1:3, k});
            % Neg_shuf = vertcat(Neg_allDay_shuf{1:3, k});
            % Pos_shuf = vertcat(Neg_allDay_shuf{1:3, k});
        elseif Day == 10
            Neg = vertcat(Neg_allDay{4:7, k});
            Pos = vertcat(Pos_allDay{4:7, k});
            % Neg_shuf = vertcat(Neg_allDay_shuf{4:7, k});
            % Pos_shuf = vertcat(Neg_allDay_shuf{4:7, k});
        end

        %
        if Day == 11 % three successful days

            if k == 1
                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];

                AA = CellListT_all_animals;
                AA(~ind_temp) = {[]};

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

                % ===== shuffle =====
                BB = CellListT_all_animals_shuf;
                BB(~ind_temp) = {[]};

                Neg_shuf = vertcat(BB{1:N_NegAnimal,:});
                Pos_shuf = vertcat(BB{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            else

                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];
                A = CellList_category_all;
                A(~ind_temp) = {[]};

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

                Neg_3trials = Neg;
                Pos_3trials = Pos;
                sgst = strcat("three successful trials");

                B = shuf_CellList_category_all;
                B(~ind_temp) = {[]};
                BB = cell(size(B));
                for i = 1:numel(B)
                    if ~isempty(B{i}) && isstruct(B{i}) && isfield(B{i}, fieldname)
                        BB{i} = B{i}.(fieldname);
                    else
                        BB{i} = [];
                    end
                end

                Neg_shuf = vertcat(BB{1:N_NegAnimal,:});
                Pos_shuf = vertcat(BB{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

            end

        end


        % === field data ===
        if isempty(Neg)
            Data_A = [];
        else
            Data_A = Neg.(fieldName);
        end
        if isempty(Pos)
            Data_B = [];
        else
            Data_B = Pos.(fieldName);
        end

        if isempty(Neg_shuf)
            Data_A_shuf = [];
        else
            Data_A_shuf = cell2mat(Neg_shuf.(fieldNameShuf));
        end
        if isempty(Pos_shuf)
            Data_B_shuf = [];
        else
            Data_B_shuf = cell2mat(Pos_shuf.(fieldNameShuf));
        end

        % z-scoring
        Data_A_shuf = (Data_A_shuf - mean(Data_A_shuf, 'omitmissing')) ./ std(Data_A_shuf, 'omitmissing');
        Data_B_shuf = (Data_B_shuf - mean(Data_B_shuf, 'omitmissing')) ./ std(Data_B_shuf, 'omitmissing');

        % === show plots ===
        fig_i = fig_indices(k);
        [p_ks(Day, k), ks2stat(Day, k), p_rank(Day, k), stats_rank{Day, k}, e_delta(Day, k)] = ...
            func_plotCdgAndBox_shuf(Data_A, Data_B, Data_A_shuf, Data_B_shuf, figrow, figcol, fig_i, fieldName, showhist);

        % sgtitle(sprintf("Day %d", Day))
        if Day<=7
            sgst = strcat("Day", num2str(Day));
        elseif Day == 8
            sgst = strcat("Day All");
        elseif Day == 9
            sgst = strcat("Day 1-3");
        elseif Day == 10
            sgst = strcat("Day 4-7");
        end
        sgtitle(sgst)
        sgst_fig{Day} = sgst;
    end


    titles = {'All cell', 'Unstable cell', 'Reliable cell', 'Non-place cell', 'Place cell'};
    for i = 1:numel(titles)
        for j = 1:2
            idx = (i-1)*2 + j;
            subplot(2,6,idx);
            ax = gca;
            ax.Title.String = {titles{i}, ax.Title.String};
        end
    end

    for i = 1:5
        subplot(2,6, i*2)
        ylabel(fieldName, 'Interpreter', 'none')
    end


end


%% For fig panel
Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

Neg = Neg_3trials;
Pos = Pos_3trials;

% === field data ===
if isempty(Neg)
    Data_A = [];
else
    Data_A = Neg.(fieldName);
end
if isempty(Pos)
    Data_B = [];
else
    Data_B = Pos.(fieldName);
end


if isempty(Neg_shuf)
    Data_A_shuf = [];
else
    Data_A_shuf = cell2mat(Neg_shuf.(fieldNameShuf));
end

if isempty(Pos_shuf)
    Data_B_shuf = [];
else
    Data_B_shuf = cell2mat(Pos_shuf.(fieldNameShuf));
end

% z-scoring
Data_A_shuf = (Data_A_shuf - mean(Data_A_shuf, 'omitmissing')) ./ std(Data_A_shuf, 'omitmissing');
Data_B_shuf = (Data_B_shuf - mean(Data_B_shuf, 'omitmissing')) ./ std(Data_B_shuf, 'omitmissing');


% close all
fig2 = figure;
% get(gcf,'Position')
set(fig2, 'Position', [1120         708         120         170])

hold on


A = Data_A;
B = Data_B;
mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

As = Data_A_shuf;
Bs = Data_B_shuf;
mAs = mean(As,'omitmissing');
mBs = mean(Bs,'omitmissing');
semAs = std(As,'omitmissing') ./ sqrt(sum(~isnan(As)));
semBs = std(Bs,'omitmissing') ./ sqrt(sum(~isnan(Bs)));

hold on
bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
ax = gca;
ax.Layer = 'top';

errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);
errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);



xticks(1:2);
xticklabels({'Neg', 'Pos'})
title('')
xlim([0.2 2.8])

disp(mA)
disp(mB)
disp(semA)
disp(semB)

sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Control: %.3f ± %.3f\n', mA, sdA);
fprintf('Casp: %.3f ± %.3f\n', mB, sdB);

% --- Wilcoxon rank-sum test (Mann–Whitney U test) ---

if isempty(Data_A) == 0 && isempty(Data_B) == 0
    [p_rank_2, ~, stats_rank_2] = ranksum(Data_A, Data_B);
    Effect = meanEffectSize(Data_A, Data_B, "Effect", "cliff" );
    e_delta = Effect.Effect;
    fprintf(sprintf('Rank-sum test p = %.3f, delta = %.3f\n', p_rank_2, e_delta));

    [p, ~, stats] = ranksum(Data_A, Data_B);
    n1 = numel(Data_A);
    n2 = numel(Data_B);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('All U1 = %.0f, U2 = %.0f, p = %.6g\n', U1, U2, p);
else

end

% try
%     [p_rank, ~, stats_rank] = ranksum(Data_A, Data_B);
%     Effect = meanEffectSize(Data_A, Data_B, "Effect", "cliff" );
%     e_delta = Effect.Effect;
%     title(sprintf('Rank-sum test p = %.3f, delta = %.3f', p_rank, e_delta));

%%

fig3 = figure;
% get(gcf,'Position')
set(fig3, 'Position', [1271         708         203         170])

hold on


hold on;
if isempty(Data_A)
    h1 = cdfplot(0);
else
    h1 = cdfplot(Data_A);
end
if isempty(Data_B)
    h2 = cdfplot(0);
else
    h2 = cdfplot(Data_B);
end
h1.Color = [0 0 0 0.8];
h1.LineWidth = 2;
h2.Color = [Ccasp 0.8];
h2.LineWidth = 2;

% h3 = cdfplot(Data_A_shuf);
if isempty(Data_A_shuf)
    h3 = cdfplot(0);
else
    h3 = cdfplot(Data_A_shuf);
end
h3.Color = [0.6 0.6 0.6 0.4];
h3.LineWidth = 2;
% h4 = cdfplot(Data_B_shuf);
if isempty(Data_B_shuf)
    h4 = cdfplot(0);
else
    h4 = cdfplot(Data_B_shuf);
end
h4.Color = [[255 200 100]/255, 0.4];
h4.LineWidth = 2;

grid off;

% --- Kolmogorov–Smirnov test ---
if isempty(Data_A) == 0 && isempty(Data_B) == 0
    try
        [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
        title(sprintf('K–S test p = %.4f', p_ks));
    catch
        p_ks = 1; ks2stat = 0;
    end
else
    p_ks = 1; ks2stat = 0;
end
title('')
ylim([0 1])
yticks([0 0.5 1])
xlabel('')
ylabel('')


end



%% func_ShowCDFandBox with shuffle data
function [fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig3] ...
    = func_ShowCDFandBoxwithShuf_2(CellListT_all_animals, CellList_category_all, CellListT_all_animals_shuf, shuf_CellList_category_all, fieldName, fieldNameShuf, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, varargin)


showhist = 0;
useField2 = false;

if nargin >= 11 && ~isempty(varargin{1})
    showhist = varargin{1};
end


%%
fig = cell(11,1);
Neg_allDay = cell(7,5);
Pos_allDay = cell(7,5);
Neg_allDay_shuf = cell(7,5);
Pos_allDay_shuf = cell(7,5);

sgst_fig    = cell(11,1);
p_ks        = nan(10,5);
ks2stat     = nan(10,5);
p_rank      = nan(10,5);
stats_rank  = cell(10,5);
e_delta     = nan(10,5);

for Day = 1:11

    fig{Day} = figure('Position', [-1584, 404, 1486, 440]);
    figrow = 2;
    figcol = 6;

    % define fields
    fields = { ...
        'CellListT_all_animals', 'Cell_unstableRatemap', ...
        'Cell_reliableRatemap', 'Cell_reliableRatemap_NotPlaceCell', ...
        'Cell_reliableRatemap_PlaceCell'};

    % fig indices
    fig_indices = [1, 3, 5, 7, 9];

    % show plots
    for k = 1:numel(fields)
        Neg_shuf = [];
        Pos_shuf = [];

        if Day<=7
            fieldname = fields{k};
            if k == 1
                A = CellListT_all_animals(:, Day);
                AA = A;

                AA_shuf = CellListT_all_animals_shuf(:, Day);
            else
                A = CellList_category_all(:, Day);

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end
            end

            % grouping Negative and Positive
            Neg = vertcat(AA{1:N_NegAnimal,:});
            Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            Neg_allDay{Day,k} = Neg;
            Pos_allDay{Day,k} = Pos;

            Neg_shuf = vertcat(AA_shuf{1:N_NegAnimal,:});
            Pos_shuf = vertcat(AA_shuf{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            Neg_allDay_shuf{Day,k} = Neg_shuf;
            Pos_allDay_shuf{Day,k} = Pos_shuf;
        end

        if Day == 8
            Neg = vertcat(Neg_allDay{:, k});
            Pos = vertcat(Pos_allDay{:, k});
        elseif Day == 9
            Neg = vertcat(Neg_allDay{1:3, k});
            Pos = vertcat(Pos_allDay{1:3, k});
        elseif Day == 10
            Neg = vertcat(Neg_allDay{4:7, k});
            Pos = vertcat(Pos_allDay{4:7, k});
        end

        if Day == 11 % three successful days
            if k == 1
                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];

                AA = CellListT_all_animals;
                AA(~ind_temp) = {[]};

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

                BB = CellListT_all_animals_shuf;
                BB(~ind_temp) = {[]};
                Neg_shuf = vertcat(BB{1:N_NegAnimal,:});
                Pos_shuf = vertcat(BB{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            else
                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];
                A = CellList_category_all;
                A(~ind_temp) = {[]};

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

                Neg_3trials = Neg;
                Pos_3trials = Pos;
                sgst = strcat("three successful trials");

                B = shuf_CellList_category_all;
                B(~ind_temp) = {[]};
                BB = cell(size(B));
                for i = 1:numel(B)
                    if ~isempty(B{i}) && isstruct(B{i}) && isfield(B{i}, fieldname)
                        BB{i} = B{i}.(fieldname);
                    else
                        BB{i} = [];
                    end
                end

                Neg_shuf = vertcat(BB{1:N_NegAnimal,:});
                Pos_shuf = vertcat(BB{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            end
        end

        % === field data ===
        if isempty(Neg)
            Data_A = [];
        else
            Data_A = Neg.(fieldName);
        end
        if isempty(Pos)
            Data_B = [];
        else
            Data_B = Pos.(fieldName);
        end

        if isempty(Neg_shuf)
            Data_A_shuf = [];
        else
            Data_A_shuf = cell2mat(Neg_shuf.(fieldNameShuf));
        end
        if isempty(Pos_shuf)
            Data_B_shuf = [];
        else
            Data_B_shuf = cell2mat(Pos_shuf.(fieldNameShuf));
        end


        % === show plots ===
        fig_i = fig_indices(k);
        [p_ks(Day, k), ks2stat(Day, k), p_rank(Day, k), stats_rank{Day, k}, e_delta(Day, k)] = ...
            func_plotCdgAndBox_shuf(Data_A, Data_B, Data_A_shuf, Data_B_shuf, figrow, figcol, fig_i, fieldName, showhist);

        if Day<=7
            sgst = strcat("Day", num2str(Day));
        elseif Day == 8
            sgst = strcat("Day All");
        elseif Day == 9
            sgst = strcat("Day 1-3");
        elseif Day == 10
            sgst = strcat("Day 4-7");
        end
        sgtitle(sgst)
        sgst_fig{Day} = sgst;
    end


    titles = {'All cell', 'Unstable cell', 'Reliable cell', 'Non-place cell', 'Place cell'};
    for i = 1:numel(titles)
        for j = 1:2
            idx = (i-1)*2 + j;
            subplot(2,6,idx);
            ax = gca;
            ax.Title.String = {titles{i}, ax.Title.String};
        end
    end

    for i = 1:5
        subplot(2,6, i*2)
        ylabel(fieldName, 'Interpreter', 'none')
    end


end

%%

%% For fig panel
Ccont = [1 1 1];
Ccasp = [255 168 60]/255;

Neg = Neg_3trials;
Pos = Pos_3trials;

% === field data ===
Data_A = Neg.(fieldName);
Data_B = Pos.(fieldName);

% Data_A_shuf = cell2mat(Neg_shuf.(fieldNameShuf));
if isempty(Neg_shuf)
    Data_A_shuf = [];
else
    Data_A_shuf = cell2mat(Neg_shuf.(fieldNameShuf));
end
% Data_B_shuf = cell2mat(Pos_shuf.(fieldNameShuf));
if isempty(Pos_shuf)
    Data_B_shuf = [];
else
    Data_B_shuf = cell2mat(Pos_shuf.(fieldNameShuf));
end


% close all
fig2 = figure;
% get(gcf,'Position')
set(fig2, 'Position', [1120         708         120         170])

hold on


A = Data_A;
B = Data_B;
mA = mean(A,'omitmissing');
mB = mean(B,'omitmissing');
semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

As = Data_A_shuf;
Bs = Data_B_shuf;
mAs = mean(As,'omitmissing');
mBs = mean(Bs,'omitmissing');
semAs = std(As,'omitmissing') ./ sqrt(sum(~isnan(As)));
semBs = std(Bs,'omitmissing') ./ sqrt(sum(~isnan(Bs)));

hold on
bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
ax = gca;
ax.Layer = 'top';

errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);
errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
    'MarkerFaceColor',Ccont, 'CapSize',5);

xticks(1:2);
xticklabels({'Neg', 'Pos'})


title('')
xlim([0.2 2.8])


% --- Wilcoxon rank-sum test (Mann–Whitney U test) ---
if isempty(Data_A) == 0 && isempty(Data_B) == 0
    [p_rank_2, ~, stats_rank_2] = ranksum(Data_A, Data_B);
    Effect = meanEffectSize(Data_A, Data_B, "Effect", "cliff" );
    e_delta = Effect.Effect;
    fprintf(sprintf('Rank-sum test p = %.3f, delta = %.3f\n', p_rank_2, e_delta));

    [p, ~, stats] = ranksum(Data_A, Data_B);
    n1 = numel(Data_A);
    n2 = numel(Data_B);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('All U1 = %.0f, U2 = %.0f, p = %.6g\n', U1, U2, p);
else
    % p_rank = 1; stats_rank = [];
    % e_delta = 0;
end

sdA = std(A,'omitmissing');
sdB = std(B,'omitmissing');
fprintf('Control: %.2f ± %.2f\n', mA, sdA);
fprintf('Casp: %.2f ± %.2f\n', mB, sdB);

%%

fig3 = figure;
% get(gcf,'Position')
set(fig3, 'Position', [1271         708         203         170])

hold on;
if isempty(Data_A)
    h1 = cdfplot(0);
else
    h1 = cdfplot(Data_A);
end
if isempty(Data_B)
    h2 = cdfplot(0);
else
    h2 = cdfplot(Data_B);
end
h1.Color = [0 0 0 0.8];
h1.LineWidth = 2;
h2.Color = [Ccasp 0.8];
h2.LineWidth = 2;

if isempty(Data_A_shuf)
    h3 = cdfplot(0);
else
    h3 = cdfplot(Data_A_shuf);
end
h3.Color = [0.6 0.6 0.6 0.4];
h3.LineWidth = 2;
if isempty(Data_B_shuf)
    h4 = cdfplot(0);
else
    h4 = cdfplot(Data_B_shuf);
end
h4.Color = [[255 200 100]/255, 0.4];
h4.LineWidth = 2;

grid off;

% --- Kolmogorov–Smirnov test ---
if isempty(Data_A) == 0 && isempty(Data_B) == 0
    try
        [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
        title(sprintf('K–S test p = %.4f', p_ks));
    catch
        p_ks = 1; ks2stat = 0;
    end
else
    p_ks = 1; ks2stat = 0;
end
title('')
ylim([0 1])
yticks([0 0.5 1])
xlabel('')
ylabel('')

end


function [p_ks, ks2stat, p_rank, stats_rank, e_delta] = func_plotCdgAndBox_shuf(Data_A, Data_B, Data_A_shuf, Data_B_shuf, figrow, figcol, fig_i, fieldName, showhist)
%% cdf plot
% showhist = 0;
% showhist = 1;

if showhist == 0
    subplot(figrow, figcol, fig_i);
    hold on;
    if isempty(Data_A)
        h1 = cdfplot(0);
    else
        h1 = cdfplot(Data_A);
    end
    if isempty(Data_B)
        h2 = cdfplot(0);
    else
        h2 = cdfplot(Data_B);
    end
    h1.Color = 'b';
    h1.LineWidth = 1;
    h2.Color = 'r';
    h2.LineWidth = 1;

    if isempty(Data_A_shuf)
        h3 = cdfplot(0);
    else
        h3 = cdfplot(Data_A_shuf);
    end
    h3.Color = 'k';
    h3.LineWidth = 1;
    if isempty(Data_B_shuf)
        h4 = cdfplot(0);
    else
        h4 = cdfplot(Data_B_shuf);
    end
    h4.Color = [1 1 1]*0.6;
    h4.LineWidth = 1;

    lgd = legend({'Negative', 'Positive', 'Neg_shuf', 'Pos_shuf'},'Location','southeast', 'Interpreter','none');
    lgd.ItemTokenSize = [10, 8];

    xlabel(fieldName, 'Interpreter','none')
    ylabel('Cumulative probability');
    grid off;

    % --- Kolmogorov–Smirnov test ---
    if isempty(Data_A) == 0 && isempty(Data_B) == 0
        try
            [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
            title(sprintf('K–S test p = %.4f', p_ks));
        catch
            p_ks = 1; ks2stat = 0;
        end
    else
        p_ks = 1; ks2stat = 0;
    end

    func_showstar_forCDF(p_ks);

elseif showhist == 1


end

%% box plot
subplot(figrow, figcol, fig_i+1)
hold on

b = boxplot([Data_A_shuf; Data_B_shuf; Data_A; Data_B], ...
    [repmat({'Neg_shuf'}, numel(Data_A_shuf), 1); ...
    repmat({'Pos_shuf'}, numel(Data_B_shuf), 1); ...
    repmat({'Negative'}, numel(Data_A), 1); ...
    repmat({'Positive'}, numel(Data_B), 1)], ...
    'Colors', [0 0 0; 0.6, 0.6, 0.6; 0 0 1; 1 0 0], 'Symbol', '', 'Whisker', 1.5);

ylabel(fieldName, 'Interpreter','none')
title('Spatial Information Comparison');
set(gca, 'Box', 'off', 'FontSize', 8);
xlim([0.2 4.8])

ax = gca;
Ymax = prctile([Data_A(:); Data_B(:)], 98);
if Ymax > ax.YLim(1)
    ax.YLim = [ax.YLim(1), Ymax*1.4];
end


try
    [p_rank, ~, stats_rank] = ranksum(Data_A, Data_B);
    Effect = meanEffectSize(Data_A, Data_B, "Effect", "cliff" );
    e_delta = Effect.Effect;
    title(sprintf('Rank-sum test p = %.3f, delta = %.3f', p_rank, e_delta));
catch
    p_rank = 1; stats_rank = [];
    e_delta = 0;
end


x = [3 4];
func_showstar_box_withShuf(p_rank, e_delta, x)

end

%
function func_showstar_box_withShuf(p_val, e_delta, x)
ax = gca;
yLim = ax.YLim;
yStar = yLim(2) + range(yLim)*0.1;
ylim([yLim(1), yStar + range(yLim)*0.3]);

p = p_val;
if p < 0.05 && abs(e_delta) > 0.147
    plot(x, [yStar yStar], 'k-', 'LineWidth', 1);

    if abs(e_delta) > 0.474
        text(mean(x), yStar + range(yLim)*0.05, '###', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    elseif abs(e_delta) <= 0.474 && abs(e_delta) > 0.33
        text(mean(x), yStar + range(yLim)*0.05, '##', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    elseif abs(e_delta) <= 0.33 && abs(e_delta) > 0.147
        text(mean(x), yStar + range(yLim)*0.05, '#', ...
            'FontSize', 12, 'HorizontalAlignment', 'center');
    end


elseif p < 0.05 && abs(e_delta) <= 0.147
    text(mean(x), yStar + range(yLim)*0.20, {'ignorable'}, ...
        'FontSize', 10, 'HorizontalAlignment', 'center');

elseif p>=0.05
    text(mean(x), yStar + range(yLim)*0.20, {'n.s.'}, ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
end

end


%% func_ShowCDFandBox with additional condition
function [fig, sgst_fig, p_ks, ks2stat, p_rank, stats_rank, e_delta, fig2, fig2_2, fig3, fig3_2] = ...
    func_ShowCDFandBox_withF2(CellListT_all_animals, CellList_category_all, FieldStats_table_animal, fieldName, fieldName2, fieldThr, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN, varargin)

showhist = 0;
% useField2 = false;

if nargin >= 11 && ~isempty(varargin{1})
    showhist = varargin{1};
end


%%
fig = cell(11,2);
sgst_fig    = cell(10,1);
p_ks        = nan(10,5);
ks2stat     = nan(10,5);
p_rank      = nan(10,5);
stats_rank  = cell(10,5);
e_delta     = nan(10,5);

Neg_allDay = cell(7,5);
Pos_allDay = cell(7,5);
Neg_FieldStats_allDay = cell(7,5);
Pos_FieldStats_allDay = cell(7,5);

F_bor_Neg_allDay = cell(7,5);
F_bor_Pos_allDay = cell(7,5);
F_cen_Neg_allDay = cell(7,5);
F_cen_Pos_allDay = cell(7,5);

%
fields = { ...
    'CellListT_all_animals', ...
    'Cell_unstableRatemap', ...
    'Cell_reliableRatemap', ...
    'Cell_reliableRatemap_NotPlaceCell', ...
    'Cell_reliableRatemap_PlaceCell'};

fig_indices = [1, 3, 5, 7, 9];
titles = {'All cell', 'Unstable cell', 'Reliable cell', 'Non-place cell', 'Place cell'};

%
for Day = 1:7
    for k = 1:numel(fields)
        categoryname = fields{k};

        if k == 1
            A = CellListT_all_animals(:, Day);
            AA = A;
        else
            A = CellList_category_all(:, Day);
            AA = cell(size(A));
            for i = 1:numel(A)
                if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, categoryname)
                    AA{i} = A{i}.(categoryname);
                else
                    AA{i} = [];
                end
            end
        end

        % --- Negative/Positive grouping ---
        Neg_allDay{Day,k} = vertcat(AA{1:N_NegAnimal,:});
        Pos_allDay{Day,k} = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

        %%
        F_border_temp = cell(N_NegAnimal+N_PosAnimal, 1);
        F_center_temp = cell(N_NegAnimal+N_PosAnimal,1);
        for animal_ind = 1:N_NegAnimal+N_PosAnimal
            A = vertcat(AA{1:N_NegAnimal+N_PosAnimal,:});
            if isempty(A) == 0
                cond_animal = A.Animal == animal_ind;
                CellInd = A.CellInd(cond_animal);
            end

            F = FieldStats_table_animal{animal_ind,:};
            if isempty(F) == 0
                F.NumPeak(F.NumPeak == 0) = 1;

                cond_cell = ismember(F.CellInd, CellInd);
                cond_area = F.Area >= 200;

                F_thrshed = F(cond_cell & cond_area,:);

                % cond_border = F_thrshed.BorderScore > 0.5;
                cond_border = F_thrshed.(fieldName2) > fieldThr;
                F_border_temp{animal_ind, 1} = F_thrshed(cond_border, :);
                F_center_temp{animal_ind, 1} = F_thrshed(~cond_border, :);
            else
                F_border_temp{animal_ind, 1} = [];
                F_center_temp{animal_ind, 1} = [];
            end
        end

        F_bor_Neg_allDay{Day,k} = vertcat(F_border_temp{1:N_NegAnimal,:});
        F_bor_Pos_allDay{Day,k} = vertcat(F_border_temp{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

        F_cen_Neg_allDay{Day,k} = vertcat(F_center_temp{1:N_NegAnimal,:});
        F_cen_Pos_allDay{Day,k} = vertcat(F_center_temp{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

    end
end

%% three successful trials
F_bor_Neg_3days = cell(numel(fields), 1);
F_bor_Pos_3days = cell(numel(fields), 1);
F_cen_Neg_3days = cell(numel(fields), 1);
F_cen_Pos_3days = cell(numel(fields), 1);

for k = 1:numel(fields)
    categoryname = fields{k};

    if k == 1
        ind_temp = [ind_NegNonNaN;ind_PosNonNaN];

        AA = CellListT_all_animals;
        AA(~ind_temp) = {[]};

        Neg = vertcat(AA{1:N_NegAnimal,:});
        Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
    else

        ind_temp = [ind_NegNonNaN;ind_PosNonNaN];
        A = CellList_category_all;
        A(~ind_temp) = {[]};


        AA = cell(size(A));
        for i = 1:numel(A)
            if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, categoryname)
                AA{i} = A{i}.(categoryname);
            else
                AA{i} = [];
            end
        end

        Neg = vertcat(AA{1:N_NegAnimal,:});
        Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

        Neg_3trials = Neg;
        Pos_3trials = Pos;
        sgst = strcat("three successful trials");

    end

    %%
    F_border_temp = cell(N_NegAnimal+N_PosAnimal, 1);
    F_center_temp = cell(N_NegAnimal+N_PosAnimal,1);

    F3 = FieldStats_table_animal;
    F3(~ind_temp) = {[]};

    for animal_ind = 1:N_NegAnimal+N_PosAnimal
        A = vertcat(AA{1:N_NegAnimal+N_PosAnimal,:});
        cond_animal = A.Animal == animal_ind;
        CellInd = A.CellInd(cond_animal);

        rowCells = F3(animal_ind,:);
        rowCells = rowCells(~cellfun(@isempty,rowCells));

        if isempty(rowCells)
            F = table();
        else
            F = vertcat(rowCells{:});
        end

        if isempty(F) == 0
            F.NumPeak(F.NumPeak == 0) = 1;

            cond_cell = ismember(F.CellInd, CellInd);
            cond_area = F.Area >= 200;

            F_thrshed = F(cond_cell & cond_area,:);

            cond_border = F_thrshed.(fieldName2) > fieldThr;
            F_border_temp{animal_ind, 1} = F_thrshed(cond_border, :);
            F_center_temp{animal_ind, 1} = F_thrshed(~cond_border, :);
        else
            F_border_temp{animal_ind, 1} = [];
            F_center_temp{animal_ind, 1} = [];
        end
    end

    F_bor_Neg_3days{k} = vertcat(F_border_temp{1:N_NegAnimal,:});
    F_bor_Pos_3days{k} = vertcat(F_border_temp{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

    F_cen_Neg_3days{k} = vertcat(F_center_temp{1:N_NegAnimal,:});
    F_cen_Pos_3days{k} = vertcat(F_center_temp{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

end


%%
close all

for ShowBorder = 1:2
    %%
    for Day = 1:11
        % fig{Day, ShowBorder} = figure('Position', [-1584, 404, 1486, 440], 'Visible','off');
        fig{Day, ShowBorder} = figure('Position', [-1584, 404, 1486, 440], 'Visible','on');
        figrow = 2; figcol = 6;

        for k = 1:numel(fields)

            if Day <= 7
                Neg_bor = F_bor_Neg_allDay{Day,k};
                Pos_bor = F_bor_Pos_allDay{Day,k};
                Neg_cen = F_cen_Neg_allDay{Day,k};
                Pos_cen = F_cen_Pos_allDay{Day,k};
            elseif Day == 8
                Neg_bor = vertcat(F_bor_Neg_allDay{:,k});
                Pos_bor = vertcat(F_bor_Pos_allDay{:,k});
                Neg_cen = vertcat(F_cen_Neg_allDay{:,k});
                Pos_cen = vertcat(F_cen_Pos_allDay{:,k});
            elseif Day == 9
                Neg_bor = vertcat(F_bor_Neg_allDay{1:3,k});
                Pos_bor = vertcat(F_bor_Pos_allDay{1:3,k});
                Neg_cen = vertcat(F_cen_Neg_allDay{1:3,k});
                Pos_cen = vertcat(F_cen_Pos_allDay{1:3,k});
            elseif Day == 10
                Neg_bor = vertcat(F_bor_Neg_allDay{4:7,k});
                Pos_bor = vertcat(F_bor_Pos_allDay{4:7,k});
                Neg_cen = vertcat(F_cen_Neg_allDay{4:7,k});
                Pos_cen = vertcat(F_cen_Pos_allDay{4:7,k});
            end

            if Day == 11 %% three successful trials
                Neg_bor = F_bor_Neg_3days{k};
                Pos_bor = F_bor_Pos_3days{k};
                Neg_cen = F_cen_Neg_3days{k};
                Pos_cen = F_cen_Pos_3days{k};
            end


            %%
            Data_A_bor = Neg_bor.(fieldName);
            Data_B_bor = Pos_bor.(fieldName);
            Data_A_cen = Neg_cen.(fieldName);
            Data_B_cen = Pos_cen.(fieldName);

            fig_i = fig_indices(k);

            if ShowBorder == 1
                [~, ~, ~, ~, ~] = ...
                    func_plotCdgAndBox_F2(Data_A_bor, Data_A_cen, Data_B_bor, Data_B_cen, figrow, figcol, fig_i, fieldName, showhist);

            elseif ShowBorder ==2
                Data_A = [Data_A_bor; Data_A_cen];
                Data_B = [Data_B_bor; Data_B_cen];
                [~, ~, ~, ~, ~] = ...
                    func_plotCdgAndBox(Data_A, Data_B, figrow, figcol, fig_i, fieldName, showhist);
            end
        end

        % === title ===
        switch Day
            case {1,2,3,4,5,6,7}
                sgst = sprintf('Day %d', Day);
            case 8
                sgst = 'All Days (1–7)';
            case 9
                sgst = 'Days 1–3';
            case 10
                sgst = 'Days 4–7';
            case 11
                sgst = 'successful three days';
        end
        sgtitle(sgst);
        sgst_fig{Day} = sgst;

        for i = 1:numel(titles)
            for j = 1:2
                idx = (i-1)*2 + j;
                subplot(2,6,idx);
                % ax = gca;
                % ax.Title.String = {titles{i}, ax.Title.String};
            end
        end

        for i = 1:5
            subplot(2,6, i*2)
            ylabel(fieldName, 'Interpreter', 'none')
        end
    end

    %%

    %%
    %% For fig panel
    Ccont = [1 1 1];
    Ccasp = [255 168 60]/255;

    k = 5; % place cell

    % three successful trials
    Neg_bor = F_bor_Neg_3days{k};
    Pos_bor = F_bor_Pos_3days{k};
    Neg_cen = F_cen_Neg_3days{k};
    Pos_cen = F_cen_Pos_3days{k};


    % === field data ===
    Data_A_bor = Neg_bor.(fieldName);
    Data_B_bor = Pos_bor.(fieldName);
    Data_A_cen = Neg_cen.(fieldName);
    Data_B_cen = Pos_cen.(fieldName);

    Data_A = [Data_A_bor; Data_A_cen];
    Data_B = [Data_B_bor; Data_B_cen];


    % close all
    fig2 = figure;
    % get(gcf,'Position')
    set(fig2, 'Position', [1120         708         120         170])

    hold on


    A = Data_A;
    B = Data_B;
    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
    semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

    hold on
    bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
    bar(2, mB, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
    ax = gca;
    ax.Layer = 'top';

    errorbar(1, mA, semA, 'Color','k', 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);
    errorbar(2, mB, semB, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);

    xticks(1:2);
    xticklabels({'Neg', 'Pos'})

    title('')
    xlim([0.2 2.8])

    fprintf('all cell\n')
    disp(mA)
    disp(mB)
    disp(semA)
    disp(semB)

    sdA = std(A,'omitmissing');
    sdB = std(B,'omitmissing');
    fprintf('Control: %.3f ± %.3f\n', mA, sdA);
    fprintf('Casp: %.3f ± %.3f\n', mB, sdB);


    try
        [p_rank, ~, stats_rank] = ranksum(Data_A, Data_B);
        Effect = meanEffectSize(Data_A, Data_B, "Effect", "cliff" );
        e_delta = Effect.Effect;

        n1 = numel(Data_A);
        n2 = numel(Data_B);
        [p, ~, stats] = ranksum(Data_A, Data_B);
        % Mann-Whitney U
        U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
        U2 = n1*n2 - U1;                            % Group B
        fprintf('All U1 = %.0f, U2 = %.0f, p = %.6g, n1 = %d, n2 = %d\n', U1, U2, p, n1, n2);
        stats_rank.U1 = U1;
        stats_rank.U2 = U2;
    catch
        p_rank = 1; stats_rank = [];
        e_delta = 0;
    end

    %%

    % close all
    fig2_2 = figure;
    % get(gcf,'Position')
    set(fig2_2, 'Position', [453   698   216   170])

    hold on

    Ccasp_b = [255 200 100]/255;

    A = Data_A_bor;
    B = Data_B_bor;
    mA = mean(A,'omitmissing');
    mB = mean(B,'omitmissing');
    semA = std(A,'omitmissing') ./ sqrt(sum(~isnan(A)));
    semB = std(B,'omitmissing') ./ sqrt(sum(~isnan(B)));

    sdA = std(A,'omitmissing');
    sdB = std(B,'omitmissing');
    fprintf('Control: %.3f ± %.3f\n', mA, sdA);
    fprintf('Casp: %.3f ± %.3f\n', mB, sdB);


    Ac = Data_A_cen;
    Bc = Data_B_cen;
    mAc = mean(Ac,'omitmissing');
    mBc = mean(Bc,'omitmissing');
    semAc = std(Ac,'omitmissing') ./ sqrt(sum(~isnan(Ac)));
    semBc = std(Bc,'omitmissing') ./ sqrt(sum(~isnan(Bc)));

    sdA = std(Ac,'omitmissing');
    sdB = std(Bc,'omitmissing');
    fprintf('Control: %.3f ± %.3f\n', mAc, sdA);
    fprintf('Casp: %.3f ± %.3f\n', mBc, sdB);


    hold on
    bar(3.5, mAc, 0.5, 'FaceColor', Ccont, 'EdgeColor', 'k');
    bar(4.5, mBc, 0.5, 'FaceColor', Ccasp, 'EdgeColor', 'none');
    bar(1, mA, 0.5, 'FaceColor', Ccont, 'EdgeColor', [0.6 0.6 0.6]);
    bar(2, mB, 0.5, 'FaceColor', Ccasp_b, 'EdgeColor', 'none');
    ax = gca;
    ax.Layer = 'top';

    errorbar(3.5, mAc, semAc, 'Color','k', 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);
    errorbar(4.5, mBc, semBc, 'Color',Ccasp*0.8, 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);
    errorbar(1, mA, semA, 'Color', [0.6 0.6 0.6], 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);
    errorbar(2, mB, semB, 'Color',Ccasp_b*0.8, 'LineWidth',0.5, ...
        'MarkerFaceColor',Ccont, 'CapSize',5);

    xticks([1 2 3.5 4.5]);
    xticklabels({'Nb', 'Pb', 'Nc', 'Pc'})
    title('')

    fprintf('center/boundary \n')
    disp(mAc)
    disp(mBc)
    disp(semAc)
    disp(semBc)

    disp(mA)
    disp(mB)
    disp(semA)
    disp(semB)

    fprintf('boundary \n')
    [p_rank, ~, stats_rank_bor] = ranksum(Data_A_bor, Data_B_bor);
    % stats_rank_bor
    % p_rank

    Effect = meanEffectSize(Data_A_bor, Data_B_bor, "Effect", "cliff" );
    e_delta = Effect.Effect

    %
    n1 = numel(Data_A_bor);
    n2 = numel(Data_B_bor);
    [p, ~, stats] = ranksum(Data_A_bor, Data_B_bor);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('border U1 = %.0f, U2 = %.0f, p = %.6g\n', U1, U2, p);
    stats_rank.U1 = U1;
    stats_rank.U2 = U2;


    [p, ~, stats] = ranksum(Data_A_bor, Data_B_bor);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('All U1 = %.0f, U2 = %.0f, p = %.6g, n1 = %d, n2 = %d\n', U1, U2, p, n1, n2);
    stats_rank.U1 = U1;
    stats_rank.U2 = U2;


    %%
    fprintf('center \n')
    [p_rank, ~, stats_rank_cen] = ranksum(Data_A_cen, Data_B_cen);
    % stats_rank_cen
    % p_rank

    Effect = meanEffectSize(Data_A_cen, Data_B_cen, "Effect", "cliff" );
    e_delta = Effect.Effect;

    %
    n1 = numel(Data_A_cen);
    n2 = numel(Data_B_cen);
    [p, ~, stats] = ranksum(Data_A_cen, Data_B_cen);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('center U1 = %.0f, U2 = %.0f, p = %.6g\n', U1, U2, p);
    stats_rank.U1 = U1;
    stats_rank.U2 = U2;

    [p, ~, stats] = ranksum(Data_A_cen, Data_B_cen);
    % Mann-Whitney U
    U1 = n1*n2 + n1*(n1+1)/2 - stats.ranksum;  % Group A
    U2 = n1*n2 - U1;                            % Group B
    fprintf('All U1 = %.0f, U2 = %.0f, p = %.6g, n1 = %d, n2 = %d\n', U1, U2, p, n1, n2);
    stats_rank.U1 = U1;
    stats_rank.U2 = U2;

    %%
    fig3 = figure;
    % get(gcf,'Position')
    set(fig3, 'Position', [1271         708         203         170])

    hold on;
    if isempty(Data_A)
        h1 = cdfplot(0);
    else
        h1 = cdfplot(Data_A);
    end
    if isempty(Data_B)
        h2 = cdfplot(0);
    else
        h2 = cdfplot(Data_B);
    end
    h1.Color = [0 0 0 0.8];
    h1.LineWidth = 2;
    h2.Color = [Ccasp 0.8];
    h2.LineWidth = 2;

    grid off;

    % --- Kolmogorov–Smirnov test ---
    try
        [~,p_ks,ks2stat] = kstest2(Data_A, Data_B);
        fprintf(sprintf('K–S test p = %.4f, stats = %.4f\n', p_ks, ks2stat));
    catch
        p_ks = 1; ks2stat = 0;
    end
    title('')
    ylim([0 1])
    yticks([0 0.5 1])
    xlabel('')
    ylabel('')

    %%
    fig3_2 = figure;
    % get(gcf,'Position')
    set(fig3_2, 'Position', [1271         708         203         170])

    hold on;
    h1 = cdfplot(Data_A_cen);
    h2 = cdfplot(Data_B_cen);
    h1.Color = [0 0 0 0.8];
    h1.LineWidth = 2;
    h2.Color = [Ccasp 0.8];
    h2.LineWidth = 2;

    h3 = cdfplot(Data_A_bor);
    h3.Color = [0.6 0.6 0.6 0.4];
    h3.LineWidth = 2;
    h4 = cdfplot(Data_B_bor);
    h4.Color = [[255 200 100]/255, 0.4];
    h4.LineWidth = 2;

    grid off;

    % --- Kolmogorov–Smirnov test ---
    try
        [~,p_ks,ks2stat] = kstest2(Data_A_bor, Data_B_bor);
        fprintf(sprintf('K–S test p = %.4f, stats = %.4f\n', p_ks, ks2stat));
        [~,p_ks,ks2stat] = kstest2(Data_A_cen, Data_B_cen);
        fprintf(sprintf('K–S test p = %.4f, stats = %.4f\n', p_ks, ks2stat));
    catch
        p_ks = 1; ks2stat = 0;
    end
    % title(sprintf('K–S test p = %.4f', p_ks));
    title('')
    ylim([0 1])
    yticks([0 0.5 1])
    xlabel('')
    ylabel('')

end


end


function [p_ks, ks2stat, p_rank, stats_rank, e_delta] = func_plotCdgAndBox_F2(Data_A_bor, Data_A_cen, Data_B_bor, Data_B_cen, figrow, figcol, fig_i, fieldName, showhist)
%% cdf plot
% showhist = 0;
% showhist = 1;

if isempty(Data_A_bor)
    Data_A_bor = 0;
end
if isempty(Data_A_cen)
    Data_A_cen = 0;
end
if isempty(Data_B_bor)
    Data_B_bor = 0;
end
if isempty(Data_B_cen)
    Data_B_cen = 0;
end

if showhist == 0
    subplot(figrow, figcol, fig_i);
    hold on;
    h1 = cdfplot(Data_A_bor);
    h2 = cdfplot(Data_B_bor);
    h1.Color = [0.5, 0.5, 1];
    h1.LineWidth = 1;
    h2.Color = [1, 0.5, 0.5];
    h2.LineWidth = 1;

    h3 = cdfplot(Data_A_cen);
    h4 = cdfplot(Data_B_cen);
    h3.Color = [0 0 1];
    h3.LineWidth = 1;
    h4.Color = [1 0 0];
    h4.LineWidth = 1;


    lgd = legend({'Neg_border', 'Pos_border', 'Neg_center', 'Pos_center'},'Location','southeast', 'Interpreter','none');
    lgd.ItemTokenSize = [10, 8];

    xlabel(fieldName, 'Interpreter','none')
    ylabel('Cumulative probability');
    grid off;

    % --- Kolmogorov–Smirnov test ---
    [~,p_ks,ks2stat] = kstest2(Data_A_cen, Data_B_cen);
    title(sprintf('K–S test p = %.4f', p_ks));

    func_showstar_forCDF(p_ks);

elseif showhist == 1


end

%% box plot
% figure;
subplot(figrow, figcol, fig_i+1)
hold on

b = boxplot([Data_A_bor; Data_B_bor; Data_A_cen; Data_B_cen], ...
    [repmat({'Neg_bor'}, numel(Data_A_bor), 1); ...
    repmat({'Pos_bor'}, numel(Data_B_bor), 1); ...
    repmat({'Neg_cen'}, numel(Data_A_cen), 1); ...
    repmat({'Pos_cen'}, numel(Data_B_cen), 1)], ...
    'Colors', [0.5 0.5 1; 1 0.5 0.5; 0 0 1; 1 0 0], 'Symbol', '', 'Whisker', 1.5);

% ylabel('Spatial information');
ylabel(fieldName, 'Interpreter','none')
title('Spatial Information Comparison');
set(gca, 'Box', 'off', 'FontSize', 8);
xlim([0.2 4.8])

ax = gca;
Ymax = prctile([Data_A_bor(:); Data_B_bor(:); Data_A_cen(:); Data_B_cen(:)], 98);
ax.YLim = [ax.YLim(1), Ymax*1.4];


% --- Wilcoxon rank-sum test (Mann–Whitney U test) ---
[p_rank, ~, stats_rank] = ranksum(Data_A_cen, Data_B_cen);
Effect = meanEffectSize(Data_A_cen, Data_B_cen, "Effect", "cliff" );
e_delta = Effect.Effect;

x = [3 4];
func_showstar_box_withShuf(p_rank, e_delta, x)


[p_rank_b, ~, stats_rank] = ranksum(Data_A_bor, Data_B_bor);
Effect = meanEffectSize(Data_A_bor, Data_B_bor, "Effect", "cliff" );
e_delta_b = Effect.Effect;
x = [1 2];
func_showstar_box_withShuf(p_rank_b, e_delta_b, x)

title({sprintf('Border Rank-sum test p = %.3f, delta = %.3f', p_rank_b, e_delta_b),...
    sprintf('Center Rank-sum test p = %.3f, delta = %.3f', p_rank, e_delta)});
end



%% Show speed mod percentage

function [fig, sgst_fig, pval] = func_CellPercentage(CellListT_all_animals, CellList_category_all, fieldName, fieldThr, Ymax, N_NegAnimal, N_PosAnimal, ind_NegNonNaN, ind_PosNonNaN)
%%

RateNeg_Day = nan(10,5);
RatePos_Day = nan(10,5);
countA = cell(10,5);
countB = cell(10,2);

for Day = 1:11

    % define fields
    fields = { ...
        'CellListT_all_animals', 'Cell_unstableRatemap', ...
        'Cell_reliableRatemap', 'Cell_reliableRatemap_NotPlaceCell', ...
        'Cell_reliableRatemap_PlaceCell'};

    varNames = {'All cells', 'unstable cells', 'Reliable rate map cell', 'non place cell', 'place cell'};
    % show plots
    for k = 1:numel(fields)
        if Day<=7
            fieldname = fields{k};
            if k == 1
                A = CellListT_all_animals(:, Day);
                AA = A;
            else
                A = CellList_category_all(:, Day);

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end
            end

            % grouping Negative and Positive
            Neg = vertcat(AA{1:N_NegAnimal,:});
            Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

            Neg_allDay{Day,k} = Neg;
            Pos_allDay{Day,k} = Pos;

        end

        if Day == 8
            Neg = vertcat(Neg_allDay{:, k});
            Pos = vertcat(Pos_allDay{:, k});
        elseif Day == 9
            Neg = vertcat(Neg_allDay{1:3, k});
            Pos = vertcat(Pos_allDay{1:3, k});
        elseif Day == 10
            Neg = vertcat(Neg_allDay{4:7, k});
            Pos = vertcat(Pos_allDay{4:7, k});
        end


        %%
        if Day == 11 % three successful days

            if k == 1
                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];

                AA = CellListT_all_animals;
                AA(~ind_temp) = {[]};

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});
            else
                ind_temp = [ind_NegNonNaN;ind_PosNonNaN];
                A = CellList_category_all;
                A(~ind_temp) = {[]};

                AA = cell(size(A));
                for i = 1:numel(A)
                    if ~isempty(A{i}) && isstruct(A{i}) && isfield(A{i}, fieldname)
                        AA{i} = A{i}.(fieldname);
                    else
                        AA{i} = [];
                    end
                end

                Neg = vertcat(AA{1:N_NegAnimal,:});
                Pos = vertcat(AA{N_NegAnimal+1 : N_NegAnimal + N_PosAnimal,:});

                Neg_3trials = Neg;
                Pos_3trials = Pos;
                sgst = strcat("three successful trials");
            end

        end
        %%


        try
            % === field data ===
            condA = Neg.(fieldName) < fieldThr;
            condB = Pos.(fieldName) < fieldThr;

            countA{Day, k} = [sum(condA), sum(~condA)];
            countB{Day, k} = [sum(condB), sum(~condB)];

            RateNeg = mean(condA) * 100;
            RatePos = mean(condB) * 100;

            RateNeg_Day(Day, k) = RateNeg;
            RatePos_Day(Day, k) = RatePos;
        catch
            RateNeg_Day(Day, k) = NaN;
            RatePos_Day(Day, k) = NaN;

        end
    end
end

%%
% === show plots ===
% close all
fig = figure('Position', [-1584, 70, 846, 774]);
pval = nan(5, 10);  % 7 days + All

for k = 1:5
    % ==== Day 1–7 ====
    subplot(5,7, (k-1)*7 + [1 3]);
    hold on
    b = bar([RateNeg_Day(1:7,k), RatePos_Day(1:7,k)]);
    set(b(1), 'FaceColor', [0.5 0.5 1]);
    set(b(2), 'FaceColor', [1 0.5 0.5]);
    title(sprintf('%s, %s', varNames{k}, fieldName), 'Interpreter','none');
    xlabel('Day'); ylabel('Ratio'); ylim([0 Ymax]); box off

    % Chi² test + star plot
    for Day = 1:7
        pval(k,Day) = chi2test(countA{Day,k}, countB{Day,k});
        addPstars_SpeedCellRatio(pval(k,Day), Day, ylim);
    end

    % ==== Day All ====
    subplot(5,7,(k-1)*7 + 4); hold on
    Day = 8;
    b1 = bar(1, RateNeg_Day(Day,k), 0.6, 'FaceColor', [0.5 0.5 1]);
    b2 = bar(2, RatePos_Day(Day,k), 0.6, 'FaceColor', [1 0.5 0.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Neg','Pos'});
    xlim([0.2 2.8]); ylim([0 Ymax]);
    ylabel('Ratio'); box off;
    title('Day all (1-7)', 'Interpreter','none');

    % Chi² test + star plot
    pval(k,Day) = chi2test(countA{Day,k}, countB{Day,k});
    addPstars_SpeedCellRatio(pval(k,Day), 1.5, ylim);


    % ==== Day 1-3 ====
    subplot(5,7,(k-1)*7 + 5); hold on
    Day = 9;
    b1 = bar(1, RateNeg_Day(Day,k), 0.6, 'FaceColor', [0.5 0.5 1]);
    b2 = bar(2, RatePos_Day(Day,k), 0.6, 'FaceColor', [1 0.5 0.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Neg','Pos'});
    xlim([0.2 2.8]); ylim([0 Ymax]);
    ylabel('Ratio'); box off;
    title('Day 1-3', 'Interpreter','none');

    % Chi² test + star plot
    pval(k,Day) = chi2test(countA{Day,k}, countB{Day,k});
    addPstars_SpeedCellRatio(pval(k,Day), 1.5, ylim);


    % ==== Day 4-7 ====
    subplot(5,7,(k-1)*7 + 6); hold on
    Day = 10;
    b1 = bar(1, RateNeg_Day(Day,k), 0.6, 'FaceColor', [0.5 0.5 1]);
    b2 = bar(2, RatePos_Day(Day,k), 0.6, 'FaceColor', [1 0.5 0.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Neg','Pos'});
    xlim([0.2 2.8]); ylim([0 Ymax]);
    ylabel('Ratio'); box off;
    title('Day 4-7', 'Interpreter','none');

    % Chi² test + star plot
    pval(k,Day) = chi2test(countA{Day,k}, countB{Day,k});
    addPstars_SpeedCellRatio(pval(k,Day), 1.5, ylim);


    % ==== successful three days ====
    subplot(5,7,(k-1)*7 + 7); hold on
    Day = 11;
    b1 = bar(1, RateNeg_Day(Day,k), 0.6, 'FaceColor', [0.5 0.5 1]);
    b2 = bar(2, RatePos_Day(Day,k), 0.6, 'FaceColor', [1 0.5 0.5]);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'Neg','Pos'});
    xlim([0.2 2.8]); ylim([0 Ymax]);
    ylabel('Ratio'); box off;
    title('successful 3 days', 'Interpreter','none');

    % Chi² test + star plot
    pval(k,Day) = chi2test(countA{Day,k}, countB{Day,k});
    addPstars_SpeedCellRatio(pval(k,Day), 1.5, ylim);
end

sgst_fig = sprintf('Speed Cell Ratio  %s', fieldName);
sgtitle(sgst_fig, 'Interpreter','none');

end


% helper functions
function p = chi2test(A, B)
tbl = [A; B];
rowSum = sum(tbl,2);
colSum = sum(tbl,1);
n = sum(tbl,'all');
expected = rowSum * colSum / n;
chi2 = sum((tbl - expected).^2 ./ expected,'all');
df = prod(size(tbl)-1);
p = 1 - chi2cdf(chi2, df);
end

function addPstars_SpeedCellRatio(p, xpos, yLim)
yStar = yLim(2)*0.9;
if p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = 'n.s.';
end
text(xpos, yStar + range(yLim)*0.05, stars, ...
    'FontSize', 8, 'HorizontalAlignment', 'center');
end




%% extract three sucessive trials

function [out, mask] = first3nonNaN(M)

[nRow, nCol] = size(M);

out = NaN(nRow, 3);
mask = false(nRow, nCol);

for i = 1:nRow
    validIdx = find(~isnan(M(i,:)));
    n = min(3, numel(validIdx));

    if n > 0
        out(i,1:n) = M(i, validIdx(1:n));
        mask(i, validIdx(1:n)) = true;
    end
end

end