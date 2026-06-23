close all; clear all

addpath(pwd);


%% Enter dir for each animal
s=1;
t=1;
group1_Dir{s,t}='path\for\group1_animal1\trial1';t=t+1;
group1_Dir{s,t}='path\for\group1_animal1\trial2';t=t+1;
group1_Dir{s,t}='path\for\group1_animal1\trial3';

s=2;
t=1;
group1_Dir{s,t}='path\for\group1_animal2\trial1';t=t+1;
group1_Dir{s,t}='path\for\group1_animal2\trial2';t=t+1;
group1_Dir{s,t}='path\for\group1_animal2\trial3';
:



s=1;
t=1;
group2_Dir{s,t}='path\for\group2_animal1\trial1';t=t+1;
group2_Dir{s,t}='path\for\group2_animal1\trial2';t=t+1;
group2_Dir{s,t}='path\for\group2_animal1\trial3';

s=2;
t=1;
group2_Dir{2,1}='path\for\group2_animal2\trial1';t=t+1;
group2_Dir{2,2}='path\for\group2_animal2\trial2';t=t+1;
group2_Dir{2,3}='path\for\group2_animal2\trial3';
:


filename = 'Data.mat';

if exist(filename, 'file')
    save('Data.mat','-append','group1_Dir','group2_Dir');
else
    save('Data.mat','group1_Dir','group2_Dir');
end

