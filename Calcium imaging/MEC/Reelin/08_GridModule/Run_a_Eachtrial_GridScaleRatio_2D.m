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

%%
s = 1; gp1 = [];
while isempty(SampleName{1,s}) == 0
    t = 1;
    while isempty(GROUP1{s,t}) == 0 && t < 4
        tmp = GROUP1{s,t}(GROUP1{s,t}(:,3) > 0, 1:13);
        gp1{s,t} = tmp;
        t = t + 1;
    end
    s = s + 1;
end

s = 1; gp2 = [];
while isempty(SampleName{2,s}) == 0
    t = 1;
    while isempty(GROUP2{s,t}) == 0 && t < 4
        tmp = GROUP2{s,t}(GROUP2{s,t}(:,3) > 0, 1:13);
        gp2{s,t} = tmp;
        t = t + 1;
    end
    s = s + 1;
end

%% Define Grid ratio based on peaks on 2D Scale vs Orientation
gp1_GSratio = nan(5,3); gp1_GSpeak = [];
for s = 1:5
    for t = 1:3
        [gp1_GSratio(s,t), gp1_GSpeak{s,t}] = Fn_2D_gSclOri(CDir, gp1{s,t}, SampleName{1,s}, s, t);
    end
end

gp2_GSratio = nan(7,3); gp2_GSpeak = [];
for s = 1:7
    for t = 1:3
        [gp2_GSratio(s,t), gp2_GSpeak{s,t}] = Fn_2D_gSclOri(CDir, gp2{s,t}, SampleName{2,s}, s, t);
    end
end

save("GridMod.mat", 'gp1_GSratio', 'gp1_GSpeak', 'gp2_GSratio', 'gp2_GSpeak')
