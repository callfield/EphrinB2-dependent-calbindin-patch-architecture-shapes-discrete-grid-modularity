
close all; clear all

load('../Data.mat');

addpath(pwd)


CDir=pwd;
load(strcat(CDir,"/Spatial_Info.mat"))
% Spatial_Info_WT & Spatial_Info_EB2
% 1: confitional entoropy (bit/sec)
% 2: confitional entoropy (bit/spike)
% 3: Normalised info

mkdir allcell
%% 
   

All_Info_WT=cell(7,4);;
All_Info_EB2=cell(7,4);;
% AllInfo=[ZgScore, Bscore, HD_Score_Bin6(:,1), SpatialInfo(:,3)];
wt=[];
for s=1:5
    for t=1:3
        DIR=wt_Dir{s,t};
        SpatialInfo=Spatial_Info_WT{s,t};
        [All_Info_WT{s,t}]=Visualise_spatial_info(CDir,DIR,SpatialInfo);
    end
end

eb=[];
for s=1:7
    for t=1:3
        if s==1&t==3
        else
            DIR=eb_Dir{s,t};
            SpatialInfo=Spatial_Info_EB2{s,t};
            Visualise_spatial_info(CDir,DIR,SpatialInfo);
            [All_Info_EB2{s,t}]=Visualise_spatial_info(CDir,DIR,SpatialInfo);
        end
    end
end
save("All_Info.mat","All_Info_EB2","All_Info_WT")

% AllInfo=[ZgScore, Bscore, HD_Score_Bin6(:,1), SpatialInfo(:,3)];
wt=[];
for s=1:5
    for t=1:3

            tmp=All_Info_WT{s,t};
            wt=[wt;tmp];
        
    end
end

eb=[];
for s=1:7
    for t=1:3

            tmp=All_Info_EB2{s,t};
            eb=[eb;tmp];
        
    end
end
writematrix(wt(:,4),"NormSpatial_Info.xlsx","Sheet","WT")
writematrix(eb(:,4),"NormSpatial_Info.xlsx","Sheet","EB2")

length(find(wt(:,1)<2& wt(:,2)<0.5 & wt(:,3)<0.1 & wt(:,4)>10))

wt=[];
for s=1:5
    for t=1:3

            tmp=All_Info_WT{s,t};
            wt=[wt;tmp];
        
    end
    length(find(wt(:,1)<2& wt(:,2)<0.5 & wt(:,3)<0.1 & wt(:,4)>10))/length(wt)

end