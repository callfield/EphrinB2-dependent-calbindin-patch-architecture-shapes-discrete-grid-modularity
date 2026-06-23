close all;
clear
% addpath(pwd)


%% analysis start
% mice_str = {'p13', 'p14','p15', 'p16'};
% mice_str_Folder = {'P13', 'P14','P15', 'P16'};

mice_str = {'n12', 'n14'};
mice_str_Folder = {'N12', 'N14'};

%% parameter that have used in pre-calculation

binSize_cm = 5;
factor = 5;

%%
for mice_ind = 1:2 
    %%
    baseFolder = 'K:\experiments\260121 Ca imaging\Behavior';

    vidFolder = fullfile(baseFolder, mice_str_Folder{mice_ind});

    subdirs = dir(vidFolder);
    A = {subdirs.name}';
    A = A( ~ismember(A, {'.', '..'}) );  % Remove . and ..


    mp4_files = A( contains(A, '.mp4') );
    mkv_files = A( contains(A, '.mkv') );
    avi_files = A( contains(A, '.avi') );

    vidfile_list = [avi_files; mp4_files;  mkv_files];


    %% load original arena ROI

    for vid_i = 1:length(vidfile_list)
        close all

        vidname = vidfile_list{vid_i};

        Ind = regexp(vidname, '[dD]ay', 'once');
        DayName = vidname(Ind:Ind + 3);

        DayName = ['Day', DayName(4)];
 

        %% load animal tracing
        s_savedfolder = 'K:\experiments\260121 Ca imaging\Behavior\TrackResults_260204';
        
        s = strcat('TrackPos_CalB', mice_str{mice_ind}, '_OF', DayName, '.mat');
        AnimalTrack = load(fullfile(s_savedfolder, s));

        T = AnimalTrack.all_trk;
        Time = T(:,1)'; Pos_x = T(:,2); Pos_y = T(:,3);
        EpochMeanV = AnimalTrack.meanVreshape'; Ts = Time(1); [val,idx]=min(abs(Time - Ts)); FrameInd = (idx:length(Time)+idx-1);
        FrameInd_aligned = (0:length(Time)-1); dt = diff(Time); Time_aligned = 0:dt(1): dt(1) * (length(Time)-1); %alined t=0 at Ca imaging start time
        s = {'frame_ind', 'original video time', 'time_aligned', 'x', 'y', 'Velocity', 'V_1secMean'};
        Table_video = [FrameInd', Time', Time_aligned', Pos_x, Pos_y, AnimalTrack.all_trk(:,4), EpochMeanV];
        Table_video = array2table(Table_video, 'VariableNames', s);

        cropRect = AnimalTrack.TrackData.cropRect;
        psh = AnimalTrack.TrackData.psh;
        psw = AnimalTrack.TrackData.psw;
        xDim = cropRect(3)*psw;
        yDim = cropRect(4)*psh;

        % Define x and y bin edges based on fixed dimensions and desired bin size.
        xEdges = 0:binSize_cm:ceil(xDim/ binSize_cm)*binSize_cm;
        yEdges = 0:binSize_cm:ceil(yDim / binSize_cm)*binSize_cm;


        %% show trace
        mousePos        = [Table_video.x, Table_video.y];

        mousePos_inFineMap = mousePos/binSize_cm * factor;

        x = mousePos_inFineMap(:,1);
        y = mousePos_inFineMap(:,2);


        %% load video

        disp(['Processing: ', vidname])

        % Load video
        vidObj = VideoReader( fullfile(vidFolder, vidname) );
        interval = 60*5; %5min
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


        %%

        close all

        fig = figure;
        imshow(MIP, []);
        if contains(vidname, '.avi')
            clim([0 0.2])
        end
        axis on

        % figure
        hold on
        x_inVid = x/psw + cropRect(1);
        y_inVid = y/psh + cropRect(2);
        plot(x_inVid, y_inVid)

        validIdx = ~isnan(x_inVid) & ~isnan(y_inVid);
        x_valid = x_inVid(validIdx);
        y_valid = y_inVid(validIdx);

        if numel(x_valid) >= 3
            k = convhull(x_valid, y_valid);
            hullPts = [x_valid(k), y_valid(k)];
        else
            hullPts = []; 
        end

        hold on
        plot(hullPts(:,1), hullPts(:,2))

        rectangle('Position', cropRect, 'EdgeColor','w')

        %%
        roi = drawpolygon;
        RoiPos = roi.Position;

        %%
        RoiPos_enclose = [RoiPos; RoiPos(1,:)];

        fig2 = figure;
        imshow(MIP, []);
        if contains(vidname, '.avi')
            clim([0 0.2])
        end
        axis on
        hold on
        plot(x_inVid, y_inVid)
        plot(hullPts(:,1), hullPts(:,2))
        rectangle('Position', cropRect, 'EdgeColor','w', 'LineStyle', '--')

        plot(RoiPos_enclose(:,1), RoiPos_enclose(:,2), 'w-', 'LineWidth', 1)


        %%
        outputDir = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\Control\redefine arena border\RedefinedROIs';

        [folder, name, ext] = fileparts(vidname);
        exportgraphics(fig2, fullfile(outputDir, strcat("Arena_Redefined_", name, ".jpg")));

        save(fullfile(outputDir, strcat("Arena_Redefined_", name, ".mat")), 'cropRect', "RoiPos_enclose", "psh", "psw");



    end
end



