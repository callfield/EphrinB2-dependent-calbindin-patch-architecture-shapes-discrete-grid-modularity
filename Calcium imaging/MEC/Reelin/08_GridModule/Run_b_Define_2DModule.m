clear all;

addpath('function')

load("GridMod.mat")
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

%% gp1_GSpeak, gp2_GSpeak
gp1_Mod{5,3} = []; gp1_ModScale{5,3} = []; gp1_ModOri{5,3} = [];
for s = 1:5
    for t = 1:3
        if size(gp1{s,t},1) > 10
            [gp1_Mod{s,t}, gp1_ModScale{s,t}, gp1_ModOri{s,t}] = Fn_Define_2DModule(CDir, ...
                gp1{s,t}, SampleName{1,s}, s, t, gp1_GSpeak{s,t});
        end
    end
end

gp2_Mod = []; gp2_ModScale = []; gp2_ModOri = [];
for s = 1:7
    for t = 1:3
        if size(gp2{s,t},1) > 10
            [gp2_Mod{s,t}, gp2_ModScale{s,t}, gp2_ModOri{s,t}] = Fn_Define_2DModule(CDir, ...
                gp2{s,t}, SampleName{2,s}, s, t, gp2_GSpeak{s,t});
        end
    end
end

% Remove one-cell modules.
for s = 1:5
    for t = 1:3
        if ~isempty(gp1_Mod{s,t})
            if min(size(gp1_Mod{s,t}{1},1), size(gp1_Mod{s,t}{2},1)) < 2
                gp1_Mod{s,t} = [];
                gp1_ModScale{s,t} = [];
                gp1_ModOri{s,t} = [];
                gp1_GSratio(s,t) = NaN;
            end
        end
    end
end

for s = 1:7
    for t = 1:3
        if ~isempty(gp2_Mod{s,t})
            if min(size(gp2_Mod{s,t}{1},1), size(gp2_Mod{s,t}{2},1)) < 2
                gp2_Mod{s,t} = [];
                gp2_ModScale{s,t} = [];
                gp2_ModOri{s,t} = [];
                gp2_GSratio(s,t) = NaN;
            end
        end
    end
end

save("GridMod.mat", 'gp1_Mod', 'gp2_Mod', 'gp1_ModScale', 'gp1_GSratio', ...
    'gp2_ModScale', 'gp1_ModOri', 'gp2_ModOri', 'gp2_GSratio', ...
    '-append')
