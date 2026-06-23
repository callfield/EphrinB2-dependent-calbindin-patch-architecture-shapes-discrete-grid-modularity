close all;
clear
% addpath(pwd)

%% make video list
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




%% Visualise trace during rec


for i = 1:length(vidname_full_paths)
    vidname = vidname_full_paths(i);


    [folder, name, ext] = fileparts(vidname);
    disp(name)


    A = load(fullfile(folder, "Track", strcat("cropRect_", name, ".mat")),"cropRect");
    C = load( fullfile(folder, "Light", strcat("RecFrames_2_", name, ".mat")));


    [thresh, metric] = multithresh(C.MaxIntensity(:,3), 1);
    Ind = C.MaxIntensity(:,3) > prctile(C.MaxIntensity(:,3), 90);
    MI = mean(C.MaxIntensity(Ind, 3));
    if metric > 0.8
        thr = thresh;
        I_start = find(C.MaxIntensity(:,3)>thr, 1) - 1;
        I_end = find(C.MaxIntensity(:,3)>thr, 1, 'last');
        if I_start == 0
            I_start = 1;
        end
        New_RecFrames = [I_start, I_end, I_end - I_start +1, C.MaxIntensity(I_end,2) - C.MaxIntensity(I_start,2)];
    elseif C.MaxIntensity(1, 3) + 50 < MI
        thr = C.MaxIntensity(1, 3) + 50;
        I_start = find(C.MaxIntensity(:,3)>thr, 1) - 1;
        I_end = find(C.MaxIntensity(:,3)>thr, 1, 'last');
        New_RecFrames = [I_start, I_end, I_end - I_start +1, C.MaxIntensity(I_end,2) - C.MaxIntensity(I_start,2)];
    else
        I_start = C.RecFrames(1,1);
        I_end = C.RecFrames(end,2);
        New_RecFrames = [I_start, I_end, I_end - I_start +1, C.MaxIntensity(I_end,2) - C.MaxIntensity(I_start,2)];
    end


    TrackData = load(fullfile(folder, "Track", strcat("Track_", name, ".mat")));
    raw_trk = TrackData.trk;


    psh = TrackData.psh;
    psw = TrackData.psw;

    RecFrames_refine = New_RecFrames;
    cropRect = A.cropRect;

    [all_trk, fig_v, fig_track, meanVreshape, fig_meanV] ...
        = Fn_ShowTrace_250910(psw, psh, RecFrames_refine, cropRect, raw_trk, name);


    %%
    OutFolder =  'K:\experiments\260121 Ca imaging\Behavior\TrackResults_260204';
    if ~exist(OutFolder, 'dir')
        mkdir(OutFolder);
    end

    s = strcat("Track ",name, ".jpg");
    exportgraphics(fig_track, fullfile(OutFolder, s),'Resolution',300)

    s = strcat("Velocity ",name, ".jpg");
    exportgraphics(fig_v, fullfile(OutFolder, s),'Resolution',300)
    
    s = strcat("Velocity_Epoch ",name, ".jpg");
    exportgraphics(fig_meanV, fullfile(OutFolder, s),'Resolution',300)

    outname  = regexprep(char(name), '^(.*?)[dD]ay(\d+).*$', '$1Day$2.mat');
    
    s = strcat("TrackPos_", outname);
    save(fullfile(OutFolder, s), 'all_trk', 'meanVreshape', 'TrackData', '-v7.3');
    close all

    clear('outname')
end
close all




%% functions

%% show trace
function [all_trk, fig_v, fig_track, meanVreshape, fig_meanV] = Fn_ShowTrace_250910(psw, psh, RecFrames_refine, cropRect, raw_trk, name)

%% skip large differences
% [folder, name, ext] = fileparts(vidname);
% disp(name)

% ylim([0 cropRect(4)*psh]);
% xlim([0 cropRect(3)*psw]);
% % ylim([cropRect(2), cropRect(4)]*psh);
% % xlim([cropRect(1), cropRect(3)]*psw);

trk=cell(size(RecFrames_refine,1),1);
for i=1:size(RecFrames_refine,1)

    str = RecFrames_refine(i,1);
    ed = RecFrames_refine(i,2);

    XY_temp = raw_trk(str:ed,2:3);

    % 
    xMin = 0;
    xMax = cropRect(3)*psw;
    yMin = 0;
    yMax = cropRect(4)*psh;

    outIdx = XY_temp(:,1) < xMin | XY_temp(:,1) > xMax | ...
        XY_temp(:,2) < yMin | XY_temp(:,2) > yMax;

    XY_temp(outIdx,:) = NaN;

    % ====================================

    V = raw_trk(str:ed,4);
    SpeedTh=60;

    % V = medfilt1(V, 3);

    maxIter = 100;
    iter = 0;

    while any(V > SpeedTh) && iter < maxIter
        [c, XY_fixed] = fixTooFast(XY_temp, V, SpeedTh);
        XY = XY_fixed;

        dt = 1/30;
        x = XY(:,1);
        y = XY(:,2);
        spx = diff(x)/dt; spx = [0; spx];
        spy = diff(y)/dt; spy = [0; spy];
        V = sqrt(spx.^2 + spy.^2);
        iter = iter + 1;
    end

    % ===== 
    xMin = 0;
    xMax = cropRect(3)*psw;
    yMin = 0;
    yMax = cropRect(4)*psh;

    outIdx = XY(:,1) < xMin | XY(:,1) > xMax | ...
        XY(:,2) < yMin | XY(:,2) > yMax;

    XY(outIdx,:) = NaN;
    % ====================================

    for d = 1:2
        x = XY(:,d);
        firstValid = find(~isnan(x), 1, 'first');
        nanIdx = isnan(x);   
        if nanIdx(end)
            lastValid = find(~nanIdx, 1, 'last');
        end

        if ~isempty(firstValid)
            x(firstValid:lastValid) = fillmissing(x(firstValid:lastValid), 'linear');
        end

        XY(:,d) = x;
        
    end

    % median filter
    XY(:,1) = medfilt1(XY(:,1), 15, 'omitnan', 'truncate');
    XY(:,2) = medfilt1(XY(:,2), 15, 'omitnan', 'truncate');

    x = XY(:,1);
    y = XY(:,2);
    spx = diff(x)/dt; spx = [0; spx];
    spy = diff(y)/dt; spy = [0; spy];
    V = sqrt(spx.^2 + spy.^2);

    trk{i}(:,1) = raw_trk(str:ed,1);
    trk{i}(:,2:3) = XY;
    trk{i}(:,4) = V;
end
% save( fullfile(folder, "RecTrack", strcat("RecTrack_", name, ".mat")), 'trk');



%% visualize
close all

fig_v = figure;
%get(gcf,'Position')
pos = [52   644   cropRect(3)   158];
set(gcf,'Position',pos);

all_trk=[];
for i=1:size(RecFrames_refine,1)
    all_trk= [all_trk;trk{i}];
end
h = plot(all_trk(:,1), all_trk(:,4), 'k');
xlabel('time [sec]')
ylabel('velocity [cm/sec]')
% clf
title(name, 'Interpreter', 'none')
xlim([all_trk(1,1), all_trk(end,1)])
% ylim([0 60])


fig_track = figure('Position', [50, 50, cropRect(3), cropRect(4)]);  % [x, y, width, height]
colormap(jet)
hold on
for i=1:size(RecFrames_refine,1)
    cplot(trk{i}(:,2),trk{i}(:,3),trk{i}(:,4),'LineWidth',1.5);
end


clim([0 20])
colorbar;
title(name, 'Interpreter', 'none')
% axis equal
axis image

set(gca, 'YDir','reverse')
ylim([0 cropRect(4)*psh]);
xlim([0 cropRect(3)*psw]);

%%

Vmed = medfilt1(V, 30); %1 sec

figure('Position', [89   851   522   126])
h = plot(all_trk(:,1)/60, Vmed, 'k');
xlabel('time [min]')
ylabel('velocity [cm/sec]')
xlim([all_trk(1,1), all_trk(end,1)]/60)
h.LineWidth = 1.0;


%% select 1s epochs to separate animal states

v2 = Vmed;
tq = all_trk(:,1);

epochT = 1; %15[s]
epoch = round(epochT/dt);
n_epoch = floor(length(v2)/epoch);
reshapedV = reshape(v2(1:n_epoch * epoch), [epoch, n_epoch]);
meanV = mean(reshapedV);
meanVreshape = reshape(repmat(meanV, epoch, 1), [1, n_epoch*epoch]);

tepoch = tq;
if length(meanVreshape) < length(tq)
    N = length(tq) - length(meanVreshape);
    Ave = mean(v2(end-N+1:end));
    meanVreshape = [meanVreshape, ones(1, N)*Ave];

end


fig_meanV = figure('Position', [89   351   1522   526]);
subplot(311)
h = plot(tq/60, v2, 'k');
xlabel('time [min]')
ylabel('velocity [cm/sec]')
hold on
xlim([tq(1) tq(end)]/60)
plot(tepoch/60, meanVreshape, 'r')
title(name,'interpreter','none')

subplot(312)
h = plot(tq/60, v2, 'k');
xlabel('time [min]')
ylabel('velocity [cm/sec]')
hold on
ylim([0 3])
yticks([0 0.2 1.0 2.0 3.0])
h = plot([tq(1) tq(end)]/60, [2 2], '--k');
h = plot([tq(1) tq(end)]/60, [1 1], '--k');
h = plot([tq(1) tq(end)]/60, [0.2 0.2], '--k');
xlim([tq(1) tq(end)]/60)
plot(tepoch/60, meanVreshape, 'r')

subplot(313)
h = plot(tq/60, v2, 'k');
xlabel('time [min]')
ylabel('velocity [cm/sec]')
hold on
ylim([0 0.5])
yticks([0 0.2 0.5 1.0 2.0])
h = plot([tq(1) tq(end)]/60, [0.2 0.2], '--k');
xlim([tq(1) tq(end)]/60)
plot(tepoch/60, meanVreshape, 'r')


%%

XY_temp = XY;
XY_temp(meanVreshape<1, :) = NaN;
C = trk{1}(:,4);
C(meanVreshape<1, :) = NaN;

fig_track = figure('Position', [50, 50, cropRect(3), cropRect(4)]);  % [x, y, width, height]
colormap(jet)
hold on

plot(XY(:,1),XY(:,2),'Color',[1 1 1]*0.7, 'LineWidth',1.5);
cplot(XY_temp(:,1),XY_temp(:,2),C,'LineWidth',1.5);


clim([0 20])
colorbar;
title(strcat(name, " Running"),'interpreter','none')
% axis equal
axis image
xlim([0 cropRect(3)*psw]);

set(gca, 'YDir','reverse')
ylim([0 cropRect(4)*psh]);
xlim([0 cropRect(3)*psw]);
end

%%
function [c, XY_fixed] = fixTooFast(XY, V, SpeedTh)

% Thr_V = 50;
Thr_V = 60;

XY_fixed = XY;
N = size(XY, 1);

% Identify too fast rows
% tooFastIdx = find(V > 100); 
tooFastIdx = find(V > Thr_V); 

% Initialize a logical index array
badIdx = false(size(XY_fixed,1),1);

% Mark the frames around each tooFastIdx (±1) as bad
for i = 1:length(tooFastIdx)
    idxRange = max(tooFastIdx(i)-1, 1):min(tooFastIdx(i), size(XY_fixed,1));
    badIdx(idxRange) = true;
end

% Set those rows to NaN
XY_fixed(badIdx, :) = NaN;

% Separate x and y coordinates
x = XY_fixed(:,1);
y = XY_fixed(:,2);
dt = 1/30;  % Frame interval

% Find valid (non-NaN) indices
validIdx = find(~any(isnan(XY_fixed), 2));

% Initialize a new speed array
V_recalc = NaN(size(x));
dt = 1/30;  % Frame time interval
maxIter = 100;  % Maximum allowed iterations
iter = 0;

while true
    % --- Get valid (non-NaN) indices ---
    validIdx = find(~any(isnan(XY_fixed), 2));

    % --- Extract x and y positions ---
    x = XY_fixed(:,1);
    y = XY_fixed(:,2);

  
    % --- Recalculate speed considering skipped NaNs ---
    V_recalc = NaN(size(x));  % Initialize speed array
    
    for i = 1:length(validIdx) - 1
        idx1 = validIdx(i);
        idx2 = validIdx(i+1);
    
        % Compute distance between valid points
        dx = x(idx2) - x(idx1);
        dy = y(idx2) - y(idx1);
        dist = sqrt(dx^2 + dy^2);
    
        % Use actual frame gap (including skipped NaNs) for time difference
        frame_diff = idx2 - idx1;
        time_diff = frame_diff * dt;
    
        % Assign the calculated speed to all intermediate frames
        V_recalc(idx1:idx2-1) = dist / time_diff;
    end


    % --- Find indices where speed exceeds threshold ---
    % tooFastIdx = find(V_recalc > 100);
    tooFastIdx = find(V_recalc > Thr_V);
    
    if isempty(tooFastIdx) || iter >= maxIter
        break;  % Stop if no fast points or iteration limit reached
    end

    % --- Mark frames around fast movement as invalid ---
    badIdx = false(size(XY_fixed,1),1);
    for i = 1:length(tooFastIdx)
        idxRange = max(tooFastIdx(i)-1, 1):min(tooFastIdx(i), size(XY_fixed,1));
        badIdx(idxRange) = true;
    end
    XY_fixed(badIdx,:) = NaN;

    iter = iter + 1;
end


% Identify too fast rows
tooFastIdx = find(V > SpeedTh); 

% Initialize a logical index array
badIdx = false(size(XY_fixed,1),1);

% Mark the frames around each tooFastIdx (±1) as bad
for i = 1:length(tooFastIdx)
    idxRange = max(tooFastIdx(i)-1, 1):min(tooFastIdx(i)+1, size(XY_fixed,1));
    badIdx(idxRange) = true;
end

% Set those rows to NaN
XY_fixed(badIdx, :) = NaN;

% Interpolate for each column independently
t          = (1:N).';                     % sample positions
goodIdx    = ~isnan(XY_fixed(:,1));   % rows not NaN (same for both cols)


c= [];

end
