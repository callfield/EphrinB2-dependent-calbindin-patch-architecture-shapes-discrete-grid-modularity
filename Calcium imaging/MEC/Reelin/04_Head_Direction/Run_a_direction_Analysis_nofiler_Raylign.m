clc; clear;
addpath(pwd)
addpath('function\')
baseDir="..\01_Each_animal";

% Directory path to the CSV files from_DeepLabCut
pathCSV = "path\to\CSV_files\from_DeepLabCut";
% CSV filename prefix
preFix = "DLC_Resnet50_....csv";

folderPaths = getOcFolders(baseDir);   % Get paths of folders starting with "Oc" and their subfolders; Reelin+ recording animals are labeled with Oc.
trialNum = numel(folderPaths);
for t = 1:trialNum
    DIR = folderPaths{t};

    % For each trial, save head and body angle data into the trial's analysis folder
    % as "HeadAngle" and "BodyAngle" inside "Angle.mat".

    saveAngle(DIR, pathCSV, preFix);
    

    caFr = 10;


    [HeadBodyDirection_FR{t}, Bined_HeadBodyDirection_FR{t}, HD_Score_Bin1{t}, HD_Score_Bin6{t}] = ...
        calculateHD_Score_nofilter_RL(DIR, caFr);
end
