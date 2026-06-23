close all; clear all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);
addpath(fullfile(scriptDir, 'function'));

load(fullfile(scriptDir, "..", "Data.mat"));

CDir=scriptDir;


GridPhase_GROUP1=cell(10,4);
GridPhase_GROUP2=cell(10,4);
AllP_GridPhase_GROUP1=cell(10,4);AllP_GridPhaseMean_allG_GROUP1=cell(10,4);AllP_GridPhaseMean_simGS_GROUP1=cell(10,4);
AllP_GridPhase_GROUP2=cell(10,4);AllP_GridPhaseMean_allG_GROUP2=cell(10,4);AllP_GridPhaseMean_simGS_GROUP2=cell(10,4);


%%
% #1 similar Grid Scale,
% #2 remove grid cell has phase >0.7 (could not use for phase analysis)
% #3 ref has higher Zscore of Grid score 
%% 


for s=1:7
    for t=1:3
        cd(group1_Dir{s,t})
        [AllP_GridPhase_GROUP1{s,t}, AllP_GridPhaseMean_allG_GROUP1{s,t}, AllP_GridPhaseMean_simGS_GROUP1{s,t}]=AllPair_GridPhaseAnalysis(CDir);
        [GridPhase_GROUP1{s,t}]=GridPhase_vs_distance_sem_patch(CDir);

    end
end
save(strcat(CDir,"\Grid_phase.mat"),"AllP_GridPhase_GROUP1","AllP_GridPhaseMean_allG_GROUP1","AllP_GridPhaseMean_simGS_GROUP1","-append")
save(strcat(CDir,"\GridPhase_GROUP1.mat"),"GridPhase_GROUP1")

for s=1:7
    for t=1:3

        cd(group2_Dir{s,t})
        [AllP_GridPhase_GROUP2{s,t}, AllP_GridPhaseMean_allG_GROUP2{s,t}, AllP_GridPhaseMean_simGS_GROUP2{s,t}]=AllPair_GridPhaseAnalysis(CDir);
        [GridPhase_GROUP2{s,t}]=GridPhase_vs_distance_sem_patch(CDir);


    end
end

save(strcat(CDir,"\Grid_phase.mat"),"AllP_GridPhase_GROUP2","AllP_GridPhaseMean_allG_GROUP2","AllP_GridPhaseMean_simGS_GROUP2","-append")
save(strcat(CDir,"\GridPhase_GROUP2.mat"),"GridPhase_GROUP2")


