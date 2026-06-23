close all; clear all
addpath('function')



load('Data.mat')
CDir=pwd;

mkdir fromBregma
%%
%repeat for all animal with each group(GROUP1 and GROUP2)
caFr = 10; % Frame rate for ca2+ recording

% GROUP1
SampleName='GROUP1_AnimalName'
s=1;

t=1;
CV=-80; % lens position form post hoc analysis
cd(group1_Dir{s,t})
[Slp_GS(t) Slp_GW(t) Grids{t}]=GridBregma_ny(SampleName,caFr,CV);

t=t+1;
CV=-50; % lens position form post hoc analysis
cd(group1_Dir{s,t})
[Slp_GS(t) Slp_GW(t) Grids{t}]=GridBregma_ny(SampleName,caFr,CV);

t=t+1;
CV=-31; % lens position form post hoc analysis
cd(group1_Dir{s,t})
[Slp_GS(t) Slp_GW(t) Grids{t}]=GridBregma_ny(SampleName,caFr,CV);


cd(CDir)
close

g1=[];
g1=[g1; Grids{1}; Grids{2}; Grids{3}];
for t=1:3
    GROUP1{s,t}(:,1:13)=Grids{t};
end



%% GROUP2
SampleName='GROUP2_AnimalName'
s=1;

t=1;
CV=-62; % lens position form post hoc analysis
cd(group2_Dir{s,t})
[Slp_GS(t) Slp_GW(t) Grids{t}]=GridBregma_ny(SampleName,caFr,CV);


t=t+1;
CV=-15; % lens position form post hoc analysis
cd(group2_Dir{s,t})
[Slp_GS(t) Slp_GW(t) Grids{t}]=GridBregma_ny(SampleName,caFr,CV);


t=t+1;
CV=-40; % lens position form post hoc analysis
cd(group2_Dir{s,t})
[Slp_GS(t) Slp_GW(t) Grids{t}]=GridBregma_ny(SampleName,caFr,CV);
t=t+1;

cd(CDir)
close

g1=[];
g1=[g1; Grids{1}; Grids{2}; Grids{3}];
for t=1:3
    GROUP2{s,t}(:,1:13)=Grids{t};
end



 cd(CDir)
 save('Data.mat','GROUP1','GROUP2','-append');
 
