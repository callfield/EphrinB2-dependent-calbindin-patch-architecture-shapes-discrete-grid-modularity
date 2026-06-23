clear all;close all;


addpath(pwd)
addpath("function")


load(fullfile(pwd, "..", "Data.mat"), "SampleName", "group1_Dir", "group2_Dir"); 
% SampleName{1,s} include animalname in animal s of WT;  
% SampleName{2,s} include animalname in animal s of EB2lz, 
load(fullfile(pwd, "Island.mat"), "wt_Island", "eb_Island") % manually define island cell; 
% wt_Island{s,t} is 0, 1,2,~ 0: no-island cell, 1: island 1, 2: island 2...
% eb_Island{s,t} is 0, 1  0: no-island cell, 1: island cell

CDir=pwd;

mkdir Ratio
mkdir MeanCor

% trace 30s; sliding window 10s; 2s increment (1/5 of s-window)
sWIN=10;

%% manually define wt_Island{s,t}{m}/eb_Island{s,t}{m}: cell id for animal s, trial t, island m


%%
% SigfCorr_ratio: ratio of sig. corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
% 
% Posi_SigfCorrMean: mean of sig. positive corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
%
% Nega_SigfCorrMean: mean of sig. positive corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr



SigfCorr_ratio=cell(5,5);
Posi_SigfCorrMean=cell(5,5);
Nega_SigfCorrMean=cell(5,5);

for s=1:5
    samplename=SampleName{1,s};
    for t=1:5
          disp(t)
        if isempty(wt_Island{s,t})==0

            DIR=group1_Dir{s,t};
            Island=wt_Island{s,t};

            if isempty(SigfCorr_ratio{s,t})==1

            [SigfCorr_ratio{s,t}, Posi_SigfCorrMean{s,t},...
                Nega_SigfCorrMean{s,t}]=...
                Fn_SubCorrelation_Island(CDir,DIR,samplename,s,t,...
                Island,sWIN);
            
            save(strcat(CDir,'\Data.mat'),...
                'SigfCorr_ratio','Posi_SigfCorrMean',...
                'Nega_SigfCorrMean','-append')


            end

        end
      
    end
end
cd(CDir)
save(strcat(CDir,'\Data.mat'),...
    'SigfCorr_ratio','Posi_SigfCorrMean','Nega_SigfCorrMean','-append')


% {
EB2_SigfCorr_ratio=cell(5,3);
EB2_Posi_SigfCorrMean=cell(5,3);
EB2_Nega_SigfCorrMean=cell(5,3);


for s=1:5
    samplename=SampleName{2,s};
    for t=1:3
          disp(t)
        if isempty(eb_Island{s,t})==0

            DIR=group2_Dir{s,t};
            Island=eb_Island{s,t};
            [EB2_SigfCorr_ratio{s,t}, EB2_Posi_SigfCorrMean{s,t},...
                EB2_Nega_SigfCorrMean{s,t}]=...
            Fn_SubCorrelation_Island(CDir,DIR,samplename,s,t,...
                Island,sWIN);

            save(strcat(CDir,'\Data.mat'),...
                'EB2_SigfCorr_ratio','EB2_Posi_SigfCorrMean',...
                'EB2_Nega_SigfCorrMean','-append')



        end
    end
end
cd(CDir)
save(strcat(CDir,'\Data.mat'),...
    'EB2_SigfCorr_ratio','EB2_Posi_SigfCorrMean',...
    'EB2_Nega_SigfCorrMean','-append')
%}
