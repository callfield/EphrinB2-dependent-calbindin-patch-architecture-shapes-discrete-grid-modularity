close all; clear all
addpath(pwd)

%% make video list 
% Define the base directory path 
behavDir = 'K:\experiments\260121 Ca imaging\Behavior';

% Get list of subdirectories in behavDir that start with "N" or "P"
subdirs = dir(behavDir);
subdirs = subdirs([subdirs.isdir]);  % Keep only directories
subdirs = subdirs( ~ismember({subdirs.name}, {'.', '..'}) );  % Remove . and ..

% Filter subdirectories that start with 'N' or 'P'
target_subdirs = subdirs( startsWith({subdirs.name}, {'N', 'P'}) );

% Initialize output lists
file_names = {};
vidname_full_paths = {};

% Loop through each target subdirectory
for i = 1:length(target_subdirs)
    subdir_path = fullfile(behavDir, target_subdirs(i).name);
    
    % Get all .mp4 files in this subdirectory
    % mp4_files = dir(fullfile(subdir_path, '*.mp4'));
    mp4_files = [dir(fullfile(subdir_path,'*.mp4')); dir(fullfile(subdir_path,'*.mkv'))];
    
    % Append file names and full paths to output lists
    for j = 1:length(mp4_files)
        file_names{end+1} = mp4_files(j).name;
        vidname_full_paths{end+1} = fullfile(subdir_path, mp4_files(j).name);
    end
end

% Convert to cell arrays or string arrays if needed
vidname_full_paths = vidname_full_paths(:);  % Column cell array
vidname_full_paths = string(vidname_full_paths);%string array


%%

for i = 1:length(vidname_full_paths)
    vidname = vidname_full_paths(i);
    [cropRect, fig] = Fn_newCamera_DefineLight_large_25909(vidname);


    [folder, name, ext] = fileparts(vidname);
    disp(['Processing: ', name])
    outputDir = fullfile(folder, "Light");
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end


    % Save cropRect
    save(fullfile(outputDir, strcat("Light_cropRect_", name, ".mat")), 'cropRect');
    exportgraphics(gcf, fullfile(outputDir, strcat("Light_ROI_", name, ".jpg")), 'Resolution', 300);
    close;

    disp("Light ROI saved successfully.");

end





%% functions
function [cropRect, fig] = Fn_newCamera_DefineLight_large_25909(vidname)
% GUI to define light ROI from video by clicking center point

[folder, name, ext] = fileparts(vidname);
disp(['Processing: ', name])


% Load video
vidObj = VideoReader(vidname);
vidObj.CurrentTime = 300;

% Use first frame to display for selection
frame = readFrame(vidObj);
img = im2double(frame(:,:,1));  % grayscale
MIP = img;  % or use median if needed

% GUI for center click
% clf;
fig = figure;
imshow(MIP, []);
title("Click to define the center of the light ROI (30×80 px)");
disp("Waiting for user click...");

[x, y] = ginput(1);  % Get single click
x = round(x);
y = round(y);

% Define 30x30 rectangle centered on click
halfSize_x = 50;
halfSize_y = 50;
cropRect = [x - halfSize_x, y - halfSize_y, halfSize_x*2, halfSize_y*2];

% Ensure cropRect is within image bounds
cropRect(1) = max(cropRect(1), 1);
cropRect(2) = max(cropRect(2), 1);
if cropRect(1) + cropRect(3) > size(MIP, 2)
    cropRect(3) = size(MIP, 2) - cropRect(1);
end
if cropRect(2) + cropRect(4) > size(MIP, 1)
    cropRect(4) = size(MIP, 1) - cropRect(2);
end

% Save cropRect
% save(fullfile(outputDir, strcat("Light_cropRect_", name, ".mat")), 'cropRect');

% Optional: show and export ROI image
X = [cropRect(1), cropRect(1)+cropRect(3), cropRect(1)+cropRect(3), cropRect(1), cropRect(1)];
Y = [cropRect(2), cropRect(2), cropRect(2)+cropRect(4), cropRect(2)+cropRect(4), cropRect(2)];
hold on;
plot(X, Y, 'w-', 'LineWidth', 2)
% exportgraphics(gcf, fullfile(outputDir, strcat("Light_ROI_", name, ".jpg")), 'Resolution', 300);
% close;
% 
% disp("Light ROI saved successfully.");

end


