clear all; close all;



addpath(pwd)
addpath("function")
load(fullfile(pwd, "..", "Data.mat"), "SampleName", "group1_Dir", "group2_Dir"); 
% SampleName{1,s} include animalname in animal s of WT;  
% SampleName{2,s} include animalname in animal s of EB2lz, 
load(fullfile(pwd, "Island.mat"), "wt_Island", "eb_Island") % manually define island cell; 
% wt_Island{s,t} is 0, 1,2,~ 0: no-island cell, 1: island 1, 2: island 2...
% eb_Island{s,t} is 0, 1  0: no-island cell, 1: island cell




mkdir dFex_230609


for s=1:5
    for t=1:5
        if isempty(wt_Island{s,t})==0
            DIR=group1_Dir{s,t};
            samplename=strcat(SampleName{1,s}," T",num2str(t));
            ISLAND=wt_Island{s,t};
            Fn_meandF_Island_ex(DIR,ISLAND,samplename, s,t);
        end
    end
end




for s=1:5
    for t=1:3
        if isempty(eb_Island{s,t})==0
            
            DIR=group2_Dir{s,t};

            samplename=strcat(SampleName{2,s}," T",num2str(t));
            ISLAND=eb_Island{s,t};
            Fn_meandF_Island_ex(DIR,ISLAND,samplename, s,t);
        end
    end
end
