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



%% make bahavior track


for mouse_i = 1:length(vidname_full_paths) %9:16 %1:8 %5:length(vidname_full_paths) %7:10 %1:3 %1:8 %9:length(vidname_full_paths) %1:8
    vidname = vidname_full_paths(mouse_i);

    if  isempty( regexp(vidname, '.mp4', 'once')) == 0 || isempty( regexp(vidname, '.mkv', 'once')) == 0
        Thresould = 100; AreaTh = 600;AreaTh2 = 1500;
        if mouse_i == 1
            Thresould = 100; AreaTh = 600;AreaTh2 = 1500;
        end
        [center, psw, psh, trk, cropRect, Circularity, Areas] = Fn_VideoTrack_newCamera_250909(vidname, Thresould, AreaTh, AreaTh2);

    elseif isempty( regexp(vidname, '.avi', 'once')) == 0

    end


    [folder, name, ext] = fileparts(vidname);
    save( fullfile(folder, "Track", strcat("Track_", name, ".mat")), ...
        'center', 'psw', 'psh', 'trk', 'cropRect', 'Circularity',...
        'Areas', '-v7.3');

    close all
end






%% functions

%% VideoTrack
function [center, psw, psh, trk, cropRect, Circularity, Areas] = Fn_VideoTrack_newCamera_250909(vidname, Thresould, AreaTh, AreaTh2)
%%
[folder, name, ext] = fileparts(vidname);
disp(name)

OutFolder = fullfile(folder, "Track");
if ~exist(OutFolder, 'dir')
    mkdir(OutFolder);
end

psw = 100/757 ;%[cm/pix]
psh = 100/757 ;%[cm/pix]

vidObj = VideoReader(vidname);
numFrames=vidObj.NumFrames;
fm = vidObj.FrameRate;
dt = 1/fm;
vidHeight = vidObj.Height;
vidWidth = vidObj.Width;


%% assume background

interval = 120;
duration = vidObj.Duration;
numFramesToRead = floor(duration / interval);
vidObj.CurrentTime = 0;
tmp_video=[];
for i = 0:numFramesToRead
    t = i * interval;
    vidObj.CurrentTime = t;
    frame = readFrame(vidObj);

    img=im2double(frame(:,:,1));% comvert to grayscale to matrix
    tmp_video(:,:,i+1)=img;%  for later averaging
end
MIP = median(tmp_video, 3); %median intensity projection

disp("background acquired")

% figure
% imshow(MIP)

%% Crop Arena position


% clf;
figure
imshow(MIP, []);
load(fullfile(folder, "Track", strcat("cropRect_", name, ".mat")),"cropRect");

% Save a preview image showing the selected area
X = [cropRect(1), cropRect(1)+cropRect(3), cropRect(1)+cropRect(3), cropRect(1), cropRect(1)];
Y = [cropRect(2), cropRect(2), cropRect(2)+cropRect(4), cropRect(2)+cropRect(4), cropRect(2)];
hold on
plot(X, Y, 'w-', 'LineWidth', 2)

exportgraphics(gcf, fullfile(folder, "Track", strcat("Background_Arena_GUI_", name, ".jpg")), "Resolution", 300);
close;

% Crop the background image
MIP = imcrop(MIP, cropRect);


%% Image complement and background subtruction

M=zeros(numFrames,1); area_array=cell(numFrames,1);biarea=cell(numFrames,1);
Circularity =zeros(numFrames,1);
susp0=[];
center = NaN(numFrames,2);
k=1;vidObj.CurrentTime = 0;

showMovie = 0;
if showMovie ==1
    figure

    frame=readFrame(vidObj);
    frame = imcrop(frame, cropRect);    % crop

    img = im2double(frame(:,:,1)); %comvert to grayscale to matrix and imvert
    img2 = MIP-img; % subtract from image
    img2(img2<0) = 0;
    img3 = mat2gray(img2)*255; % covert to matrix to grayscale image

    %made binary image with dif image and animal color
    I2bi = img3;
    I2bi = zeros(size(img3));
    mask = img3 >= Thresould & frame(:,:,1) < 40 & frame(:,:,2) < 30  & frame(:,:,3) < 40 ;
    I2bi(mask)=1; %made binary image
    hImg = imshow(I2bi);
    hold on
    hPt = plot(NaN, NaN, 'r+', 'MarkerSize',10,'LineWidth',2); 

end

SE_close = strel("square",10);
SE_open = strel("square",10);

vidObj.CurrentTime = 0;
while hasFrame(vidObj)
    frame=readFrame(vidObj);
    frame = imcrop(frame, cropRect);    % crop

    img = im2double(frame(:,:,1)); %comvert to grayscale to matrix and imvert
    img2 = MIP-img; % subtract from image
    img2(img2<0) = 0;
    img3 = mat2gray(img2)*255; % covert to matrix to grayscale image

    %made binary image with dif image and animal color
    I2bi = img3;
    I2bi = zeros(size(img3));
    mask = img3 >= Thresould & frame(:,:,1) < 40 & frame(:,:,2) < 30  & frame(:,:,3) < 40 ;
    I2bi(mask)=1; %made binary image

    img_temp = imclose(I2bi, strel("rectangle",[4 4]) );  % dilating and eroding
    img_temp = imopen(img_temp, strel("rectangle",[3 3]));

    % imshow(img_temp)
    I2bi = img_temp;

    I2bi = imclose(I2bi, SE_close);
    I2bi = imopen(I2bi, SE_open);

    %tracking maximum area's particle
    cc = bwconncomp(I2bi);%label binary object
    stats = regionprops(cc, img, 'Area', 'Centroid',"Circularity","MajorAxisLength",  "MinorAxisLength", "Perimeter", ...
        'MeanIntensity', 'MaxIntensity');


    if ~isempty(stats)
        circ = vertcat(stats.Circularity);
        validIdx = find(circ > 0.5);

        if ~isempty(validIdx)
            areas = vertcat(stats(validIdx).Area);
            areaMask = (areas > AreaTh) & (areas < AreaTh2);

            if any(areaMask)
                validIdx2 = validIdx(areaMask);
                areas2 = areas(areaMask);

                [~, relative_id] = max(areas2);
                id = validIdx2(relative_id);

                % save position of centroid and area
                M(k,1) = stats(id).Area;
                c = stats(id).Centroid;
                x = round(c(1));% sometimes there are same maxima
                y = round(c(2));% sometimes there are same maxima

                xx=x-5:x+5;yy=y-5:y+5;
                xx=xx(xx>0&xx<cropRect(3));
                yy=yy(yy>0&yy<cropRect(4));
                I2bi(yy,xx)=0.7;
                center(k,:) =c; % centroid of closest area

                Circularity(k,1) = stats(id).Circularity;
            end
        end

    end

    % writeVideo(writerObj, uint8(I2bi * 255));
    k=k+1;

    if showMovie ==1

        set(hImg, 'CData', I2bi);   
        set(hPt, 'XData', c(1), 'YData', c(2));   
        title(sprintf('Mouse Position %.1f s', vidObj.CurrentTime));
        drawnow limitrate nocallbacks  % ← fast
        % pause(0.001)
    end
end

disp("Binalized")

%%
numFrames = size(center,1); % .avi files (or similar formats) where vidObj.NumFrames is not accurat

center_fixed = fixFrozenCenter(center);

%  Deifinition of suspecious frame  rescue some frame (not too big&not too fast)
x = center_fixed(:,1);
y = center_fixed(:,2);
spx = diff(x)/dt; spx = [0; spx]; %pix/s
spy = diff(y)/dt; spy = [0; spy]; %pix/s
spy = spy*psh; %cm/s
spx = spx*psw; %cm/s
v = sqrt(spx.^2 + spy.^2);



% Mark detection errer frame

susp=[];
trk = [x*psw, y*psh, v];
trk(susp,:)=nan(length(susp),3);

tq = linspace(0, vidObj.Duration-dt, numFrames);
tq = tq.';
trk = [tq, trk];

% csvwrite(strcat(folder, "\Track\Raw_Track ",name, ".csv"), trk);


%% visualise

h = plot(trk(:,1), trk(:,4), 'k');
xlabel('time [sec]')
ylabel('velocity [cm/sec]')
clf

figure('Position', [50, 50, cropRect(3), cropRect(4)]);  % [x, y, width, height]
colormap(jet)
xx=x*psw; xx(susp)=[];yy=y*psh;yy(susp)=[];vv=v;vv(susp)=[];
cplot(xx,yy,vv,'LineWidth',1.5);
ylim([0 cropRect(4)*psh]);
xlim([0 cropRect(3)*psw]);
colorbar;
exportgraphics(gcf, strcat(folder, "\Track\Track ",name, ".jpg"),'Resolution',300)
clf


% check
tmp_M=M;
tmp_M(susp)=[];
histogram(tmp_M, 150)
xlabel('Area (pic)')
ylabel('number of frame')
% exportgraphics(gcf, strcat(folder, "\Track\Nosusp_Area_",name, ".jpg"),'Resolution',300)
% % exportgraphics(gcf, strcat(folder, "\Track\Nosusp_Area_",name, ".pdf"))
clf


Areas = M;

end



%%
function center_fixed = fixFrozenCenter(center)
% Replace frozen rows (consecutive identical rows) by interpolated values
%
% center : N×2 matrix  [x  y]
% method : (optional) interpolation method for interp1
%          'linear' (default), 'pchip', 'spline', etc.
%
% center_fixed : matrix with frozen rows replaced

% if nargin < 2
%     method = 'pchip';  % default interpolation
% end

center_fixed = center;                % copy
N            = size(center,1);

% 1) Identify rows that are identical to the previous row
isSame       = [false; all(diff(center)==0, 2)];  % logical N×1
frozenIdx    = find(isSame);                      % indices to correct
if isempty(frozenIdx)
    return                                    % nothing to fix
end

% 2) Mark frozen rows as NaN (both columns)
center_fixed(frozenIdx,:) = NaN;

% 3) Interpolate for each column independently
t          = (1:N).';                     % sample positions
goodIdx    = ~isnan(center_fixed(:,1));   % rows not NaN (same for both cols)

for col = 1:2
    center_fixed(:,col) = interp1( ...
        t(goodIdx), ...
        center_fixed(goodIdx,col), ...
        t, ...
        'pchip', ...
        'extrap');    % use extrap to handle NaNs at start/end
end

end

