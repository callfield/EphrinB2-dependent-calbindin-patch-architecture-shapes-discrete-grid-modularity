clear all;close all

addpath(pwd);

addpath('function');


load('../Data.mat');

CDir=pwd;
mkdir Visualise
%% output: Border_Score
% B_Score
caFr=10;
group1_B_Score=cell(7,3);
for s=1:7
    for t=1:3
        DIR=group1_Dir{s,t};
        cd(DIR) ;
        [ group1_B_Score{s,t}]=Border_score_v2(CDir,caFr);
    end
end
 cd(DIR)
save(strcat(CDir,"\Data_Borderv2.mat"),"group1_B_Score")



group2_B_Score=cell(7,3);
for s=1:7
    for t=1:3

        DIR=group2_Dir{s,t};
        cd(DIR) ;
        [ group2_B_Score{s,t}]=Border_score_v2(CDir,caFr);

    end
end

save(strcat(CDir,"\Data_Borderv2.mat"),"group2_B_Score","-append")

 
%%
cd(CDir)
%% border cell cdf

gru1=[];
for s=1:7
    for t=1:3
        gru1=[gru1;group1_B_Score{s,t}];
    end
end
    

gru2=[];
for s=1:7
    for t=1:3
        gru2=[gru2;group2_B_Score{s,t}];
    end
end
    
% cdf plot 
data1=gru1;
data2=gru2;
name1= "GROUP1";name2= "GROUP2";
xlbl= "Border Score v2";ylbl="Proportion";ttl= "Border Score v2";
filename= "Border Score v2.jpg";OutDir="";
CdfBarplot_bd(OutDir,data1,data2,name1, name2, xlbl,ylbl,ttl,filename)

length(find(gru1>0.5))
length(find(gru1>0.5))/length(gru1)
length(gru1)

length(find(gru2>0.5))
length(find(gru2>0.5))/length(gru2)
length(gru2)

writematrix(gru1, "BorderScore_v2.xlsx","Sheet","GROUP1")
writematrix(gru2, "BorderScore_v2.xlsx","Sheet","GROUP2")