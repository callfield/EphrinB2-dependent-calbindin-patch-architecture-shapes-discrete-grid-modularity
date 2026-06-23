clear all;close all

addpath('function');

GP1_SpeedCorr2_NewP=cell(5,3);
GP2_SpeedCorr2_NewP=cell(7,3);


addpath(pwd);

load('../Data.mat');
CDir=pwd;
mkdir Speed_vs_Hz_NewP

% cd('..\AnimalData\')

%% output: GP1_SpeedCorr2_NewP, GP2_SpeedCorr2_NewP
% GP1_SpeedCorr2_NewP{R / P/ vHZ}
% R(cell, time-shifted Correlation (Velocity vs Hz))
% P(cell, pvalue of time-shifted Correlation (Velocity vs Hz))
% vHZ(velocity, cell , shifted frame)
% Pval(cell)
% nonMoveFR
%%

caFr=10;
for s=1:5
    for t=1:3
        DIR=group1_Dir{s,t};
        cd(DIR);

        [GP1_SpeedCorr2_NewP{s,t} ]=Speed_vs_meandFCorr2_NewP(CDir,caFr);
        save(strcat(CDir,"\SpeedCorr2_NewP.mat"),"GP1_SpeedCorr2_NewP","-append")
    end
end

caFr=10;
for s=1:7
    for t=1:3
        DIR=group2_Dir{s,t};
        cd(DIR);
            [GP2_SpeedCorr2_NewP{s,t} ]=Speed_vs_meandFCorr2_NewP(CDir,caFr);
            save(strcat(CDir,"\SpeedCorr2_NewP.mat"),"GP2_SpeedCorr2_NewP","-append")

        
    end
end

