close all;
clear
% addpath(pwd)

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




%% Extract Rec frame


Thresould = 100; AreaTh = 15;AreaTh2 = 25;
for i = 5:length(vidname_full_paths) %1:length(vidname_full_paths) %1:8 %length(vidname_full_paths) %1:length(vidname_full_paths)
    vidname = vidname_full_paths(i);
    Fn_Light_Duration_250909(vidname, Thresould, AreaTh, AreaTh2);
end






%% functions

function Fn_Light_Duration_250909(vidname, Thresould, AreaTh, AreaTh2)


[folder, name, ext] = fileparts(vidname);
disp(strcat(name, " extract Rec start end Frames"))

lightFolder = fullfile(folder, "Light");
if ~exist(lightFolder, 'dir')
    mkdir(lightFolder);
end

psw = 100/540 ;%[cm/pix] 100cm/520pix
psh = 100/520 ;%[cm/pix] 100cm/520pix

fprintf(strcat("read ", name, "\n"))

vidObj = VideoReader(vidname);
% frameRate = vidObj.FrameRate;
numFrames=vidObj.NumFrames;
fm = vidObj.FrameRate;
dt = 1/fm;
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;


%% assume background

interval = 60;
duration = vidObj.Duration;
numFramesToRead = floor(duration / interval);
vidObj.CurrentTime = 0;
tmp_video=[];
for i = 0%:numFramesToRead
    t = i * interval;
    vidObj.CurrentTime = t;
    frame = readFrame(vidObj);

    img=im2double(frame(:,:,1));% comvert to grayscale to matrix

    tmp_video(:,:,i+1)=img;%  for later averaging

end
MIP = median(tmp_video, 3); %median intensity projection

disp("background acquired")

%% Crop Arena position


clf;
imshow(MIP, []);
load(fullfile(folder, "Light", strcat("Light_cropRect_", name, ".mat")),"cropRect");

% Save a preview image showing the selected area
X = [cropRect(1), cropRect(1)+cropRect(3), cropRect(1)+cropRect(3), cropRect(1), cropRect(1)];
Y = [cropRect(2), cropRect(2), cropRect(2)+cropRect(4), cropRect(2)+cropRect(4), cropRect(2)];
hold on
plot(X, Y, 'w-', 'LineWidth', 2)

% exportgraphics(gcf, fullfile(folder, "Light", strcat("Light_ROI_", name, ".jpg")));
% close;

% Crop the background image
MIP = imcrop(MIP, cropRect);


%% Image complement and background subtruction

writerObj = VideoWriter( strcat(folder, "\Light\CroppedLight_",name));

writerObj.FrameRate = fm;
% open the video writer
open(writerObj);
M=zeros(numFrames,1);
k=1;vidObj.CurrentTime = 0;
N = vidObj.NumFrames;
MaxIntensity = zeros(N,3);
while hasFrame(vidObj)
    frame=readFrame(vidObj);
    frame = imcrop(frame, cropRect);    % crop

    img = im2double(frame(:,:,1)); %comvert to grayscale to matrix and imvert
    img2 = img - MIP; % subtract from image
    img2(img2<0) = 0;
    img3 = mat2gray(img2)*255; % covert to matrix to grayscale image

    MaxIntensity(k,:) = [k, vidObj.CurrentTime, max(img(:))*255];

    %made binary image with dif image and animal color
    I2bi = img3;
    I2bi = zeros(size(img3));
    mask = img3 >= Thresould & frame(:,:,1) > 180 & frame(:,:,2) > 180  & frame(:,:,3)> 180  ;
    I2bi(mask)=1; %made binary image


    %tracking maximum area's particle
    cc = bwconncomp(I2bi);  % label binary objects
    stats = regionprops(cc, img, 'Area', 'Centroid', "Circularity", ...
        "MajorAxisLength", "MinorAxisLength", "Perimeter", ...
        'MeanIntensity', 'MaxIntensity');

    if ~isempty(stats)
        circ = vertcat(stats.Circularity);
        % validIdx = find(circ > 0.5);
        validIdx = find(circ > 0.2);

        if ~isempty(validIdx)
            areas = vertcat(stats(validIdx).Area);
            [~, relative_id] = max(areas);
            id = validIdx(relative_id);

            % check if area is between 15 and 25
            selectedArea = stats(id).Area;
            M(k,1) = selectedArea;
            center(k,:) = stats(id).Centroid;
        end
    end

    % writeVideo(writerObj, uint8(I2bi * 255));
    writeVideo(writerObj, uint8(img * 255));

    if mod(k,100)==0
        fprintf("Writing video, current time: %d\n", round(vidObj.CurrentTime *10)/10)
    end

    k=k+1;
end
close(writerObj);% close the video writer

disp("Binalized")

AreaTh = 3;

% find file showing rec start end
RecFrames=[];
[RecFrames(:,1), RecFrames(:,2)] = detectAboveThresholdSegments(M, AreaTh);
RecFrames(:,3) = RecFrames(:,2) - RecFrames(:,1) +1 ;
% RecFrames
save( fullfile(folder, "Light", strcat("RecFrames_2_", name, ".mat")), 'RecFrames', 'MaxIntensity');

end


%%

function [startFrames, endFrames] = detectAboveThresholdSegments(M, threshold)
    if nargin < 2
        threshold = 10;  % default threshold
    end

    M = M(:);  % ensure column vector

    above = M >= threshold;
    changes = diff([0; above; 0]);  % add padding to catch edges

    startFrames = find(changes == 1);  % rising edge → start of segment
    endFrames = find(changes == -1) - 1;  % falling edge → end of segment
end
