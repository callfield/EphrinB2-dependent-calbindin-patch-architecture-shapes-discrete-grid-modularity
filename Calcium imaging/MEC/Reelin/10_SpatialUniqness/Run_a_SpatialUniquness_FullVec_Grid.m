% Build 50 x 50 spatial maps from selected cell sets across 10,000 frames.
% For each base bin, compute the correlation map against all other bins and
% estimate error distances (max, mean, total) for bins above each threshold.
% Repeat this process and keep the cell set with the smallest total error.

clear all;close all

addpath(pwd);
addpath("function");


load('../Data.mat');
CDir=pwd;


%%


mkdir MaxError_GridmaxFR

% Condition x animal x trial
% cond=1: Group1, cond=2: Group2
ErrDist  = cell(2, 7, 3);

% Threshold sets (r_th, p_th)
r_list = [0.2,  0.2,  0.3,  0.3,  0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.8, 0.8];
p_list = [0.05,  0.01, 0.05,  0.01, 0.05,  0.01, 0.05,  0.01, 0.05,  0.01, 0.05,  0.01];
th_tags = {'r02p05','r02p01',  'r03p05','r03p01', 'r04p05','r04p01', 'r05p05', 'r05p01', 'r06p05', 'r06p01', 'r08p05', 'r08p01'};  % labels used as field/file suffixes

nTh = numel(r_list);
if ~( numel(r_list)==numel(p_list) && numel(r_list)==numel(th_tags) )
    error('r_list, p_list, and th_tags must have the same length');
end

caFr = 10;

%% ----- condition 1: group1 -----
cond = 1;  % index for Group1
for s = 1:5
    for t = 1:3
        DIR = group1_Dir{s, t};
        DATA = GROUP1{s, t};

        numGrid = length(find(DATA(:, 3) > 0));
        if numGrid > 9
            

            DV_pos = DATA(:, 1);
            Grid_Scale = DATA(:, 3);

            [ErrDist{cond, s, t}] = Fn_SpatialUniquness_FullVec_Grid(DIR, caFr, r_list, p_list, th_tags, DV_pos, Grid_Scale);

        end
    end
end

save(strcat(CDir, "\SpatialUniquness_Data.mat"),  "ErrDist", "-append");

%% ----- condition 2: group2 -----
cond = 2;  % index for Group2

for s = 1:7
    for t = 1:3
        DIR = group2_Dir{s, t};
        DATA = GROUP2{s, t};
        numGrid = length(find(DATA(:, 3) > 0));
        if numGrid > 9
            DV_pos = DATA(:, 1);
            Grid_Scale = DATA(:, 3);
            [ErrDist{cond, s, t}] = Fn_SpatialUniquness_FullVec_Grid(DIR, caFr, r_list, p_list, th_tags, DV_pos, Grid_Scale);
        end
    end
end

save(strcat(CDir, "\SpatialUniquness_Data.mat"),  "ErrDist", "-append");
