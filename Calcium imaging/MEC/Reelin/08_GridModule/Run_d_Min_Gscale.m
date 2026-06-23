clear all; close all;

addpath(pwd);
addpath('function')

load('../Data.mat')
load('GridMod.mat')

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

gp1 = [];
for s = 1:5
    for t = 1:3
        if isempty(gp1_Mod{s,t}) == 0
            tmp = GROUP1{s,t}(GROUP1{s,t}(:,3) > 0, 1:5);
            tmp = tmp(gp1_Mod{s,t}{1},:);
            tmp = tmp(tmp(:,1) > -2200,:);
            gp1 = [gp1; tmp];
        end
    end
end

gp2 = [];
for s = 1:7
    for t = 1:3
        if isempty(gp2_Mod{s,t}) == 0
            tmp = GROUP2{s,t}(GROUP2{s,t}(:,3) > 0, 1:5);
            tmp = tmp(gp2_Mod{s,t}{1},:);
            tmp = tmp(tmp(:,1) > -2200,:);
            gp2 = [gp2; tmp];
        end
    end
end

[h, p] = ttest2(gp1(:,3), gp2(:,3))

writematrix(gp1, '0-150_mod1_GridS.xlsx', "Sheet", "gp1")
writematrix(gp2, '0-150_mod1_GridS.xlsx', "Sheet", "gp2")
