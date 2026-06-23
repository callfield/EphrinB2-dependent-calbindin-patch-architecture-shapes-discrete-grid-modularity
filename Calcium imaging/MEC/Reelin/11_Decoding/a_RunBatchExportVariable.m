clc; clearvars; close all;
addpath("function");

% Set this to the animal data root directory.
baseDir = "path\for\AnimalData";
if ~isfolder(baseDir)
    error("baseDir does not exist: %s", baseDir);
end

folderPaths = getOcFolders(baseDir);
trialNum = numel(folderPaths);
for t = 1:trialNum
    DIR = folderPaths{t};
    ExportVariable(DIR);
end
