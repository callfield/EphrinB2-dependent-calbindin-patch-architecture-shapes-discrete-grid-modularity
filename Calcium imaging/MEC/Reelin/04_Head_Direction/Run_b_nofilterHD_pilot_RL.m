clc; clear;
addpath(pwd)

CDir=pwd;
addpath('function')
load("../Data.mat");

mkdir("nofilter_HDcell");

group1=[];
group1_HD_Score=cell(7,3);
Move_group1_HD_Score=cell(7,3);
for s=1:7
    for t=1:3
        DIR=group1_Dir{s,t};
        [group1_HD_Score{s,t}, Move_group1_HD_Score{s,t}]=HD_score_nofilter_RL(DIR,CDir,SampleName{1,s}, t);
        close all

        group1=[group1;group1_HD_Score{s,t}];
    end
end
save(strcat(CDir,"\Data_HD.mat"),"group1_HD_Score","Move_group1_HD_Score");

group2=[];
group2_HD_Score=cell(7,3);
Move_group2_HD_Score=cell(7,3);
for s=1:7
    for t=1:3

        DIR=group2_Dir{s,t};
        cd(DIR) ;
        [group2_HD_Score{s,t}, Move_group2_HD_Score{s,t}]=HD_score_nofilter_RL(DIR,CDir,SampleName{2,s}, t);
        
        close all
        group2=[group2;group2_HD_Score{s,t}];
    end
end
save(strcat(CDir,"\Data_HD.mat"),"group2_HD_Score","Move_group2_HD_Score","-append");

% cdf plot 
data1=group1(:,1);


length(find( group1(:,2)<0.05  & group1(:,4)<0.05 & group1(:,6)<40  ))
length(find( group1(:,2)<0.05  & group1(:,4)<0.05 & group1(:,6)<40  ))/length(group1)
length(group1)

length(find(group2(:,2)<0.05 & group2(:,4)<0.05 & group2(:,6)<40 ))
length(find(group2(:,2)<0.05 & group2(:,4)<0.05 & group2(:,6)<40 ))/length(group2)
length(group2)

cd(CDir)
writematrix(group1, "HDScore.xlsx","Sheet","group1")
writematrix(group2, "HDScore.xlsx","Sheet","group2")
