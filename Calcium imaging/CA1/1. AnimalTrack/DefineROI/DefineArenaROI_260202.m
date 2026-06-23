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
    [cropRect, fig] = Fn_newCamera_DefineArena_250909(vidname);


    [folder, name, ext] = fileparts(vidname);
    disp(['Processing: ', name])
    outputDir = fullfile(folder, "Track");
    if ~exist(outputDir, 'dir')
        mkdir(outputDir)
    end

    %%
    save(fullfile(outputDir, strcat("cropRect_", name, ".mat")), 'cropRect');

    % Optional: save visualization of selected area
    X = [cropRect(1), cropRect(1)+cropRect(3), cropRect(1)+cropRect(3), cropRect(1), cropRect(1)];
    Y = [cropRect(2), cropRect(2), cropRect(2)+cropRect(4), cropRect(2)+cropRect(4), cropRect(2)];
    hold on;
    plot(X, Y, 'w-', 'LineWidth', 2)
    exportgraphics(gcf, fullfile(outputDir, strcat("Arena_GUI_", name, ".jpg")));
    close;
    disp("Crop region saved successfully.");

end




%% functions

function [cropRect, fig] = Fn_newCamera_DefineArena_250909(vidname)
% GUI to define arena from video and save cropRect

[folder, name, ext] = fileparts(vidname);
disp(['Processing: ', name])


% Load video
vidObj = VideoReader(vidname);
interval = 180;
duration = vidObj.Duration;
numFramesToRead = floor(duration / interval);
vidObj.CurrentTime = 0;

% Estimate background using median projection
tmp_video = [];
for i = 0:numFramesToRead
    t = i * interval;
    vidObj.CurrentTime = t;
    frame = readFrame(vidObj);
    img = im2double(frame(:,:,1));  % grayscale
    tmp_video(:,:,i+1) = img;
end
MIP = median(tmp_video, 3);

% GUI crop
% clf;
close all
fig = figure;
imshow(MIP, []);
title("Drag to select the arena region");
disp("Waiting for user to select crop area...");
cropRect = round(getrect);  % [x, y, width, height]


% Optional: save visualization of selected area
X = [cropRect(1), cropRect(1)+cropRect(3), cropRect(1)+cropRect(3), cropRect(1), cropRect(1)];
Y = [cropRect(2), cropRect(2), cropRect(2)+cropRect(4), cropRect(2)+cropRect(4), cropRect(2)];
hold on;
plot(X, Y, 'w-', 'LineWidth', 2)
% exportgraphics(gcf, fullfile(outputDir, strcat("Arena_GUI_", name, ".jpg")));
% close;

% disp("Crop region saved successfully.");

end
