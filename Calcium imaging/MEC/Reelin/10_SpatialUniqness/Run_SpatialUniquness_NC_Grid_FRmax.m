%% Run_CorrDecode_NormCell_totaletc_230804 used random picked 1000 cell set
%% this(Run_CorrDecode_NormCell_totaletc_230804 with R2 selected cells) select top 52 cell and make random ---- cell set
clear all;close all


addpath(pwd);
addpath("function");


load('../Data.mat');
CDir=pwd;


%%

mkdir MeanError_GridmaxFR 
mkdir TotalError_GridmaxFR
mkdir MaxError_GridmaxFR

zMean_ErrDist_maxFRGrid=cell(2,7,3);
zTotal_ErrDist_maxFRGrid=cell(2,7,3);
zMax_ErrDist_maxFRGrid=cell(2,7,3);

caFr=10;
for s=1:5 
    for t=1:3
                   DIR=wt_Dir{s,t};
            % xxx=split(wt_Dir{s,t},"\");
            % DIR=strcat('..\..\AnimalData\', char(xxx(end-1)), "\",char(xxx(end)));
        numGrid=length(find(WT{s,t}(:,3)>0));
        if numGrid>9
            if s==1&&t==1
                caFr=6;
            else
                caFr=10;
            end
            [Mean_ErrDist_maxFRGrid{1,s,t}, Total_ErrDist_maxFRGrid{1,s,t}, Max_ErrDist_maxFRGrid{1,s,t}]=...
                Fn_SpatialUniquness_NC_Grid_FRmax_230817(DIR,caFr);
        end

    end
end
save(strcat(CDir,"\SpatialUniquness_Data.mat"), ...
    "Mean_ErrDist_maxFRGrid","Total_ErrDist_maxFRGrid","Max_ErrDist_maxFRGrid", "-append")

caFr=10;
for s=1:7
    for t=1:3
         if s==1&t==3
         else
         numGrid=length(find(EB2{s,t}(:,3)>0));
            if numGrid>9
                        DIR=eb_Dir{s,t};
            % xxx=split(eb_Dir{s,t},"\");
            % DIR=strcat('..\..\AnimalData\', char(xxx(end-1)), "\",char(xxx(end)));
            [Mean_ErrDist_maxFRGrid{2,s,t}, Total_ErrDist_maxFRGrid{2,s,t}, Max_ErrDist_maxFRGrid{2,s,t}]=...
                Fn_SpatialUniquness_NC_Grid_FRmax_230817(DIR,caFr);
        end
         end
    end
end

save(strcat(CDir,"\SpatialUniquness_Data.mat"), ...
    "Mean_ErrDist_maxFRGrid","Total_ErrDist_maxFRGrid","Max_ErrDist_maxFRGrid", "-append")

