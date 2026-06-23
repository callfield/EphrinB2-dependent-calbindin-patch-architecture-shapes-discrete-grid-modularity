%%
clear
close all

% get(0, 'DefaultUicontrolFontName')

s =  'MS UI Gothic';
set(0, 'DefaultAxesFontName', s);
set(0, 'DefaultTextFontName', s);

%% set parameters
binSize_cm      = 5;

Spike_thr = 6; % spike detection threshold, Spike_thr * MAD
occThresh = 2; % Remove bins stayed <2s

% gaussian filter parameters
sigma           = 1.5;
kernelSize      = [5 5];


%% analysis start
% mice_str = {'p13', 'p14','p15', 'p16'};
% mice_str_Folder = {'P13', 'P14','P15', 'P16'};
mice_str = {'n12', 'n14'};
mice_str_Folder = {'N12', 'N14'};


for mice_ind = 1:2 %1:4
    for Day = 1:6 %:6
        %% load validated Ca traces
        fprintf(strcat('Mouse', mice_str_Folder{mice_ind}, " Day%d\n"), Day);

        s_CellTraceFolder = 'H:\experiments H drive\251203 Ca imaging\code_251211\CellValidation\Control\results_Validation';
        s = strcat('CellTraces_Validated_', mice_str{mice_ind}, '_d', num2str(Day), '.mat');
        if exist(fullfile(s_CellTraceFolder, s), 'file') == 0
            fprintf(strcat("File not found: ", 'Mouse', mice_str_Folder{mice_ind}, " Day%d\n"), Day);
            continue;
        end
        Trace = load(fullfile(s_CellTraceFolder, s));

        % baseline correction
        trace_smooth_all = zeros(size(Trace.interp_traces_accepted));
        window_size = 500; % baseline drift correction
        for i = 1:size(Trace.interp_traces_accepted, 2)
            interp_trace = Trace.interp_traces_accepted(:,i);
            % baseline drift correction
            baseline_drift = movmedian(interp_trace, window_size);
            trace_detrended = interp_trace - baseline_drift;
            % smoothing
            trace_smooth = smoothdata(trace_detrended, 'movmean', 5);
            trace_smooth_all(:, i) = trace_smooth;
        end
        calciumTraces = [Trace.full_ts', trace_smooth_all];

        %% load video tracking
        s_savedfolder = 'H:\experiments H drive\260121 Ca imaging\Behavior\TrackResults_260204';
        s = strcat('TrackPos_CalB', mice_str{mice_ind}, '_OFDay', num2str(Day), '.mat');
        AnimalTrack = load(fullfile(s_savedfolder, s));

        T = AnimalTrack.all_trk;
        Time = T(:,1)'; Pos_x = T(:,2); Pos_y = T(:,3);

        EpochMeanV = AnimalTrack.meanVreshape';

        Ts = Time(1);
        [val,idx]=min(abs(Time - Ts));
        FrameInd = (idx:length(Time)+idx-1);
        FrameInd_aligned = (0:length(Time)-1);
        dt = diff(Time);
        Time_aligned = 0:dt(1): dt(1) * (length(Time)-1); %alined t=0 at Ca imaging start time

        s = {'frame_ind', 'original video time', 'time_aligned', 'x', 'y', 'Velocity', 'V_1secMean'};

        Table_video = [FrameInd', Time', Time_aligned', Pos_x, Pos_y, AnimalTrack.all_trk(:,4), EpochMeanV];
        Table_video = array2table(Table_video, 'VariableNames', s);
        %Time_aligned: alined behavior video recording to Ca imaging start = 0 sec.

        VideoSampling   = numel(Time_aligned) / Time_aligned(end);

        %% Detect Spikes and align to video timestamps
        disp("Running now: SpikeExtraction")
        % tic

        thr = Spike_thr; % [threshold] = thr * MAD

        TraceData = calciumTraces;
        [CaProcessingParm, SpkTime, Peaks] = func_SpkDetection(TraceData, thr);
        [SpkTimeAligned, SpkTimeAligned_table] = func_AlignSpikeToVideoTime(Table_video, SpkTime);
        SpkDataMapped = SpkTimeAligned;

        % toc
        % Output parameters
        SpkDetectionResults.CaProcessingParm = CaProcessingParm;
        SpkDetectionResults.SpkTime = SpkTime;
        SpkDetectionResults.Peaks = Peaks;
        SpkDetectionResults.SpkTimeAligned_table = SpkTimeAligned_table;

        %% Extract Running Frames in Behavior video

        smoothedSpeed = EpochMeanV;

        % --- Step 1: Define Thresholds and Preallocate Storage
        thresholds_speed = 0:20;  % speed thresholds in cm/sec
        nThresholds_speed = length(thresholds_speed);
        idx = cell(nThresholds_speed, 1);      % cell array to store indices for each threshold

        % --- Step 2: Find Indices Where smoothedSpeed Exceeds Each Threshold
        for i = 1:nThresholds_speed
            idx{i} = find(smoothedSpeed > thresholds_speed(i));
        end

        % --- Step 3: Save All Results in a Structure
        SpeedResults = struct();
        SpeedResults.smoothedSpeed = smoothedSpeed;  % original speed vector
        SpeedResults.Thresholds = thresholds_speed;        % thresholds used
        SpeedResults.idx = idx;               % indices for each threshold
        % SpeedResults.Timestamp = datetime('now');    % optional timestamp

        edges = 0:20;
        nBins = length(edges)-1;
        SpeedCounts = histcounts(smoothedSpeed, edges);                % number of bin appearance
        frameRate = VideoSampling;
        SpeedTimeHist   = SpeedCounts / frameRate;


        %% Rate map for single cell with speed filtering
        %% outputs preallocation
        n_Cell = numel(SpkTimeAligned_table);

        FiringRate_cell = nan(n_Cell,1);

        % rate map related variables
        rateMap_ori_sp          = cell(nThresholds_speed, n_Cell);
        rateMap_thr_sp          = cell(nThresholds_speed, n_Cell);
        occMap_ori_sp           = cell(nThresholds_speed, 1);
        occMap_thr_sp           = cell(nThresholds_speed, 1);
        spikeMap_Rawcount_sp    = cell(nThresholds_speed, n_Cell);
        Edges_sp                = cell(nThresholds_speed, n_Cell);
        Dur_Running_sp          = nan(nThresholds_speed, 1);
        areaMap_sp              = cell(nThresholds_speed, n_Cell);
        runSpkIdx_sp      = cell(nThresholds_speed, n_Cell);
        neuronPos_sp      = cell(nThresholds_speed, n_Cell);
        numEvent_sp       = nan(nThresholds_speed, n_Cell);


        XY_oldnew_cell = cell(nThresholds_speed, n_Cell);
        occMapNormalized_sp = cell(5,1);
        occMapUniform_Gaussian_sp = cell(5,1);
        rateMapUniform_sp       = cell(nThresholds_speed, n_Cell);
        rateMapUniform_Gaussian_sp     = cell(nThresholds_speed, n_Cell);


        %% output of field analysis

        FRD_FieldAnalysis_Redo_cell = cell(n_Cell,1);
        FRD_ACGanalysis_cell = cell(n_Cell,1);


        %% output of boot
        nCell = n_Cell;

        HalvesCorr_Pearson_even_cell        = nan(nCell,1);
        HalvesCorr_speed_Pearson_even_cell  = nan(nCell,1);

        HalvesCorr_Pearson_Boot_cell        = cell(nCell,1);
        CI_rateMap_cell                     = nan(nCell,2);
        CI_rateMap_90_cell                     = nan(nCell,2);
        R_boot_mean_rateMap_cell            = nan(nCell,1);
        HalvesCorr_speed_Pearson_Boot_cell  = cell(nCell,1);
        CI_speed_cell                       = nan(nCell,2);
        CI_speed_90_cell                       = nan(nCell,2);
        R_boot_mean_speed_cell              = nan(nCell,1);

        HalvesCorr_Pearson_Boot_cell_shuf           = cell(nCell,1);
        p_rank_ratemap_cell                 = nan(nCell,1);
        Cliff_delta_ratemap_cell            = nan(nCell,1);

        HalvesCorr_speed_Pearson_Boot_cell_speed    = cell(nCell,1);
        p_rank_speed_cell                   = nan(nCell,1);
        Cliff_delta_speed_cell              = nan(nCell,1);


        %% outputs preallocation (latest)

        FieldAnalysis_cell  = cell(n_Cell,1);
        SpInfo_cell         = cell(n_Cell, 1);
        Coherence_cell = cell(n_Cell, 1);

        SpInfo_shuf_cell = cell(n_Cell, 1);
        Coherence_shuf_cell = cell(n_Cell, 1);

        speedModResult_cell = cell(n_Cell, 1);

        HalveCompareField_stability_cell = cell(n_Cell, 1);
        HalveCompareSpeed_stability_cell = cell(n_Cell, 1);
        HalveCompareField_reliability_cell = cell(n_Cell, 1);
        HalveCompareSpeed_reliability_cell = cell(n_Cell, 1);

        Halves_shuf_stab_cell = cell(n_Cell, 1);
        Halves_shuf_reliab_cell = cell(n_Cell, 1);

        %%

        count = 0;
        for i_Cell = 1:n_Cell
            %%
            count = count + 1;
            overallTic = tic;
            close all

            % set constant variables
            allFrame_vid = round(Table_video.frame_ind);
            mousePos        = [Table_video.x, Table_video.y];
            spikeData       = SpkTimeAligned_table{1, i_Cell};
            spikeFrames     = spikeData.frame_ind;
            spikePos        = [spikeData.x, spikeData.y];

            % tic
            for thr_sp = 1:nThresholds_speed
                %%
                % spike position during running
                runFrames = SpeedResults.idx{thr_sp};
                runSpkIdx = ismember(spikeFrames, runFrames);
                neuronPos = spikePos(runSpkIdx, :);
                mousePos_run = mousePos(runFrames,:);

                % calculate rate map using adaptive binning
                [occMapMerged, spkMapMerged, rateMapMerged, xEdges, yEdges, nBin, Dur_Running, areaMap, rateMapMerged_ori, occMapMerged_ori] = ...
                    func_adaptiveBinRateMap(mousePos_run, neuronPos, AnimalTrack, binSize_cm, VideoSampling, occThresh, sigma, kernelSize);

                % summarize results
                rateMap_ori_sp{thr_sp, i_Cell}          = rateMapMerged_ori;
                rateMap_thr_sp{thr_sp, i_Cell}      = rateMapMerged;
                occMap_ori_sp{thr_sp, 1}                = occMapMerged_ori;
                occMap_thr_sp{thr_sp, 1}                = occMapMerged;

                spikeMap_Rawcount_sp{thr_sp, i_Cell}   = spkMapMerged;
                Edges_sp{thr_sp, i_Cell} = {xEdges, yEdges};
                Dur_Running_sp(thr_sp, 1) = Dur_Running;
                areaMap_sp{thr_sp, i_Cell}   = areaMap;

                runSpkIdx_sp{thr_sp, i_Cell}        = runSpkIdx;
                neuronPos_sp{thr_sp, i_Cell}        = neuronPos;
                numEvent_sp(thr_sp, i_Cell)       = sum(spkMapMerged(:), 'omitmissing');
            end
            % toc %0.086sec

            %% Show Ca traces abd spikes
            [fig, FiringRate] = func_showCaTraceAndSpikes_251013(AnimalTrack, Table_video, SpkTimeAligned_table, TraceData, runSpkIdx_sp, i_Cell, CaProcessingParm, SpeedResults);

            FiringRate_cell(i_Cell,1) = FiringRate;

            %% Show occupancy map


            X_old_sp = cell(5,1); Y_old_sp = cell(5,1); X_new_sp = cell(5,1); Y_new_sp = cell(5,1);

            doShow = 1;
            for thr_sp = 1:5
                [rateMapUniform, rateMapUniform_Gaussian, occMapNormalized, occMapGaussian, X_old, Y_old, X_new, Y_new, occMapUniform_Gaussian] = ...
                    func_MakeAndShowGaussianRateMap(AnimalTrack, occMap_ori_sp, Edges_sp, thr_sp, i_Cell, areaMap_sp, binSize_cm, thresholds_speed, sigma, kernelSize, rateMap_thr_sp, doShow);
                rateMapUniform_sp{thr_sp, i_Cell} = rateMapUniform;
                rateMapUniform_Gaussian_sp{thr_sp, i_Cell} = rateMapUniform_Gaussian;
                occMapNormalized_sp{thr_sp, 1} = occMapNormalized;
                occMapUniform_Gaussian_sp{thr_sp, 1} =  occMapUniform_Gaussian;

                X_old_sp{thr_sp, 1} = X_old; Y_old_sp{thr_sp, 1} = Y_old; X_new_sp{thr_sp, 1} = X_new; Y_new_sp{thr_sp, 1} = Y_new;
            end


            XY_oldnew_cell{thr_sp, i_Cell} = [X_old_sp, Y_old_sp , X_new_sp, Y_new_sp];

            %% show fine rate map

            [rateMap_Fine_sp, PeakFiringRate_sp, xs_sp, ys_sp, factor] = ...
                func_ComputeAndShowFineRateMap(AnimalTrack, Edges_sp, i_Cell, spikeMap_Rawcount_sp, rateMap_thr_sp, rateMapUniform_sp, binSize_cm, rateMapUniform_Gaussian_sp);


            mouseStr = mice_str_Folder{mice_ind};
            origID = Trace.AcceptedMetric(i_Cell,1);

            s_animal = sprintf( strcat('Mouse#', mice_str_Folder{mice_ind}, " Day%d"), Day);
            s_CellID = sprintf('Cell ID #%d (orig. #%d)', i_Cell, Trace.AcceptedMetric(i_Cell,1));
            s_figtitle = strcat(s_animal, ", ", s_CellID);
            sgtitle(strcat(s_figtitle, " whole duration"))

            %% Place Field Analysis, etc. peak firing rate and field size

            thr_sp = 3; % only focus on 2cm/s

            %MAD threshold is too low for images with many zero intensity
            [fig_field, FieldStats_ratio, FieldStats_STD] = func_FieldAnalysis(rateMap_Fine_sp, thr_sp, binSize_cm, factor);

            FieldAnalysis.PeakFiringRate    =  FieldStats_STD.PeakFR(thr_sp);
            FieldAnalysis.FieldSize_Peak    =  FieldStats_STD.FieldSize_Peak(thr_sp);
            FieldAnalysis.FieldNum          =  FieldStats_STD.FieldNum(thr_sp);
            FieldAnalysis.FieldStats_ratio  =  FieldStats_ratio;
            FieldAnalysis.FieldStats_STD    =  FieldStats_STD;

            %% Spatial Info

            occMap = occMapUniform_Gaussian_sp{thr_sp, 1} * binSize_cm^2;
            occMap = occMap * binSize_cm^2; % in each bin, change [s/cm2] to [s] (not needed)
            rateMap = rateMapUniform_Gaussian_sp{thr_sp, i_Cell};

            [meanRate, info_perspike, info_persec, spatialSparsity, spatialSelectivity] = func_SkaggsSpatialInfo_251014(rateMap, occMap);

            SpInfo.meanRate         =  meanRate;
            SpInfo.info_perspike    =  info_perspike;
            SpInfo.info_persec      =  info_persec;
            SpInfo.spatialSparsity  =  spatialSparsity;
            SpInfo.spatialSelectivity  =  spatialSelectivity;

            clear('meanRate', 'info_perspike', 'info_persec', 'spatialSparsity', 'spatialSelectivity')

            %% Coherence

            rateMap = rateMapUniform_sp{thr_sp, i_Cell};
            [Cohe_Z, Cohe_SpCorr, SmoothRateMap_forCoherence] = func_SpatialCoherence_251014(rateMap);

            Coherence.Z = Cohe_Z;
            Coherence.SpCorr = Cohe_SpCorr;


            %% shuffling

            % figure
            fprintf(strcat(num2str(Spike_thr), "MAD Running now: Compute shuffling ", s_animal, " ", s_CellID, " / ", num2str(n_Cell), "\n") );
            t1 = tic;

            [SpInfo_shuf, Coherence_shuf, shuffIdx, neuronPos_shuff] = ...
                func_ShuffleTest_adapt(thr_sp, SpeedResults, spikeFrames, mousePos, AnimalTrack, occThresh, sigma, kernelSize, VideoSampling, binSize_cm, SpInfo, Coherence, spikePos);

            info_obs        = SpInfo_shuf.info_obs;
            info_shuf       = SpInfo_shuf.info_shuf;
            info_obs_persec = SpInfo_shuf.info_obs_sec;

            SpInfo_z = SpInfo_shuf.SpInfo_z;
            SpInfo_sec_z = SpInfo_shuf.SpInfo_sec_z;

            meanRate = SpInfo.meanRate;
            spatialSelectivity = SpInfo.spatialSelectivity;

            set(0,'CurrentFigure',fig)
            subplot(5,6,12)
            h = histogram(info_shuf,30); hold on
            xline(info_obs,'r','LineWidth',2);
            xlabel('Spatial information (bits/spike)');
            ylabel('Count');
            s_test = sprintf('Shuffle test (p=%.4f)', SpInfo_shuf.pval_spinfo);
            if max([info_obs; info_shuf(:)])>0
                xlim([0 max([info_obs; info_shuf(:)])*1.2])
            end
            title({'Sp Info and Suffle test'})

            subplot(5,6,18)
            s_Info = sprintf('Spatial Info = %.2f bit/spike', info_obs);
            s_mRatio = sprintf('Mean firing ratio = %.3f Hz', meanRate);
            st = strcat("At speed > ", num2str(SpeedResults.Thresholds(thr_sp)), " cm/s, ");
            s_Spz = sprintf('z-scored Spatial Info = %.2f', SpInfo_z);

            s_SpInfoSec = sprintf('Spatial Info (per time) = %.2f bit/sec', info_obs_persec);
            s_SpInfoSecTst = sprintf('Spatial Info Shuffle: p = %.4f', SpInfo_shuf.pval_spsecinfo);
            s_SpInfoSec_z = sprintf('z-scored Spatial Info (per time) = %.2f', SpInfo_sec_z);

            s_Sparsity = sprintf('sp. selectivity = %.2f', spatialSelectivity);
            s_SparsityTest = sprintf('sp. selectivity test p = %.4f', SpInfo_shuf.pval_selectivity);
            s_Cohe = sprintf('Coherence R = %.2f', Coherence.SpCorr);

            s = {st, s_mRatio, s_Info, s_test, s_Spz, s_SpInfoSec, s_SpInfoSecTst, s_SpInfoSec_z, s_Sparsity, s_SparsityTest, s_Cohe};
            tex = text(.0,.5, s); %#ok<NASGU>
            axis off

            clear('st', 's_mRatio', 's_Area', 's_Info', 's_test', 's_Spz', 's_SpInfoSec', 's_SpInfoSecTst', 's_SpInfoSec_z', 's_Sparsity', 's_SparsityTest', 's_Cohe')

            tShuffleSP = toc(t1);
            fprintf('shuffling RatioMap: %.3f s\n', tShuffleSP);


            %% speed modulation analysis

            doShuffle = 1; nShuffle = 1000; doFitlm = 1; doFitLogistic = 1;
            minBinTime = 60;
            speedModResult = func_analyzeSpeedModulation_251014(spikeFrames, smoothedSpeed, VideoSampling, doShuffle, nShuffle, doFitlm, doFitLogistic, minBinTime);

            set(0,'CurrentFigure',fig)
            subplot(5,6,24)
            bar(speedModResult.speedCenters, speedModResult.ratePerBin, 'FaceColor',[0.3 0.6 0.8]);
            xlabel('Speed (cm/s)');
            ylabel('Firing rate (Hz)');
            title(strcat("Speed mod. ", sprintf( '(r = %.2f, p = %.3f)', speedModResult.corr_r, speedModResult.corr_p)));
            xlim([0 15])

            subplot(5,6,30)
            st = strcat("Other speed modulation metrics");
            s_Pregression = sprintf('p of regression = %.3f', speedModResult.regression_basics.p_slope);
            s_Pshuffle = sprintf('p for shuffle = %.3f', speedModResult.shuffle_p);
            s_SPinfo = sprintf('Speed Info = %.2f [bit/spike]', speedModResult.speedInfo.I_spike);
            s_SPinfo_p = sprintf('info p for shuffle = %.3f', speedModResult.speedInfo.I_spikeShuf_p);

            s = {st, s_Pregression, s_Pshuffle, s_SPinfo, s_SPinfo_p};
            tex = text(.0,.5,s);
            axis off


            %% Stability, calc rate mal divided by halves of recording time
            % close all
            speedInd_list = [1 2 3 4 5]; % speed thresholds using for calculation

            t2 = tic;
            doImShow = 1; % show figure
            doCompare = 'Stability';
            [fig2, HalveCompareField, HalveCompareSpeed] = ...
                func_HalveCompare(TraceData, allFrame_vid, Time_aligned, mousePos, spikeFrames, spikePos, SpeedResults,...
                binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, nThresholds_speed, CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, speedInd_list,...
                doImShow, doCompare);
            fig_stability = fig2;
            HalveCompareField_stability = HalveCompareField;
            HalveCompareSpeed_stability = HalveCompareSpeed;
            % toc


            doCompare = 'Reliability';
            [fig2, HalveCompareField, HalveCompareSpeed] = ...
                func_HalveCompare(TraceData, allFrame_vid, Time_aligned, mousePos, spikeFrames, spikePos, SpeedResults,...
                binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, nThresholds_speed, CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, speedInd_list,...
                doImShow, doCompare);
            fig_reliability = fig2;
            HalveCompareField_reliability = HalveCompareField;
            HalveCompareSpeed_reliability = HalveCompareSpeed;

            sttime = toc(t2);
            fprintf(sprintf('stability and reliability %3f s\n', sttime))

            %% shuffle stability


            % fprintf('shuffle stability\n')
            doCompare = 'Stability';
            t2 = tic;
            Halves_shuf_stab = func_StabilityShuffle(TraceData, allFrame_vid, Time_aligned, mousePos, ...
                spikeFrames, spikePos, SpeedResults, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, ...
                nThresholds_speed, CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, doCompare, shuffIdx);
            sttime = toc(t2);
            fprintf(sprintf('shuffle stability %3f s\n', sttime))

            doCompare = 'Reliability';
            t2 = tic;
            Halves_shuf_reliab = func_StabilityShuffle(TraceData, allFrame_vid, Time_aligned, mousePos, ...
                spikeFrames, spikePos, SpeedResults, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, ...
                nThresholds_speed, CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, doCompare, shuffIdx);
            sttime = toc(t2);
            fprintf(sprintf('shuffle reliability %3f s\n', sttime))

            %% create cell array of outputs

            FieldAnalysis_cell{i_Cell,1} = FieldAnalysis;
            SpInfo_cell{i_Cell,1} = SpInfo;
            Coherence_cell{i_Cell,1} = Coherence;

            SpInfo_shuf_cell{i_Cell,1} = SpInfo_shuf;
            Coherence_shuf_cell{i_Cell,1} = Coherence_shuf;

            speedModResult_cell{i_Cell,1} = speedModResult;

            HalveCompareField_stability_cell{i_Cell,1} = HalveCompareField_stability;
            HalveCompareSpeed_stability_cell{i_Cell,1} = HalveCompareSpeed_stability;
            HalveCompareField_reliability_cell{i_Cell,1} = HalveCompareField_reliability;
            HalveCompareSpeed_reliability_cell{i_Cell,1} = HalveCompareSpeed_reliability;

            Halves_shuf_stab_cell{i_Cell,1} = Halves_shuf_stab;
            Halves_shuf_reliab_cell{i_Cell,1} = Halves_shuf_reliab;


            %% save figures

            %set save folder
            baseFolder = 'H:\experiments H drive\260121 Ca imaging\Behavior';
            saveFolder = fullfile(baseFolder, 'PlaceFieldRawAnalysis_260204', ...
                sprintf('%dMAD', thr), mice_str_Folder{mice_ind}, ...
                sprintf('Day%d', Day));
            saveFolder_Basic = fullfile(saveFolder, 'Basics_260204');
            [status,msg,msgID] = mkdir(saveFolder_Basic);

            % common name
            animalID = sprintf('Mouse%sDay%d', mice_str_Folder{mice_ind}, Day);
            cellID   = sprintf('CellID#%d', i_Cell);


            figList = {
                'Basics_',     fig;
                'Stability_',  fig_stability;
                'Reliability_', fig_reliability
                };

            % save images
            for k = 1:size(figList,1)
                fname = sprintf('%s%s_%s.jpeg', figList{k,1}, animalID, cellID);
                saveas(figList{k,2}, fullfile(saveFolder_Basic, fname));
            end


            %% field analysis re-do
            close all

            thr_sp = 3;
            [fig_field, FieldAnalysis_Redo, ACGanalysis] = func_FieldMetrics_260204_PeakRate30(binSize_cm, factor, i_Cell, rateMap_Fine_sp, thr_sp, mice_ind, Day, occMapUniform_Gaussian_sp, AnimalTrack);
            FRD_FieldAnalysis_Redo_cell{i_Cell,1} = FieldAnalysis_Redo;
            FRD_ACGanalysis_cell{i_Cell,1} = ACGanalysis;

            % save image
            OutDir_FRD = fullfile(saveFolder, 'FieldRedo_260204_Ratio');
            [status,msg,msgID] = mkdir(OutDir_FRD);

            animalID = sprintf('Mouse%sDay%d', mice_str_Folder{mice_ind}, Day);
            cellID   = sprintf('CellID#%d', i_Cell);
            fname = sprintf('FieldRedo %s_%s.jpeg', animalID, cellID);

            exportgraphics(fig_field, fullfile(OutDir_FRD, fname));

            %% Bootstrap rate map
            close all

            t2 = tic;
            [fig_boot, HalvesCorr_Pearson_even, HalvesCorr_speed_Pearson_even, HalvesCorr_Pearson_Boot, CI_rateMap, CI_rateMap_90, R_boot_mean_rateMap, ...
                HalvesCorr_speed_Pearson_Boot, CI_speed, CI_speed_90, R_boot_mean_speed, HalvesCorr_Pearson_shuf, p_rank, Cliff_delta, HalvesCorr_speed_Pearson_shuf, p_rank_speed, Cliff_delta_speed] = ...
                func_BootstrapRatemapReliability_260204(thr_sp, SpeedResults, sigma, kernelSize, mice_ind, Day, i_Cell, AnimalTrack, spikeFrames, spikePos, binSize_cm);

            sttime = toc(t2);
            fprintf(sprintf('Bootstrap %3f s\n', sttime))

            % save variables
            HalvesCorr_Pearson_even_cell(i_Cell,1)        = HalvesCorr_Pearson_even;
            HalvesCorr_speed_Pearson_even_cell(i_Cell,1)  = HalvesCorr_speed_Pearson_even;
            HalvesCorr_Pearson_Boot_cell{i_Cell,1}        = HalvesCorr_Pearson_Boot;
            CI_rateMap_cell(i_Cell,:)                     = CI_rateMap;
            CI_rateMap_90_cell(i_Cell,:)                  = CI_rateMap_90;
            R_boot_mean_rateMap_cell(i_Cell,1)            = R_boot_mean_rateMap;
            HalvesCorr_speed_Pearson_Boot_cell{i_Cell,1}  = HalvesCorr_speed_Pearson_Boot;
            CI_speed_cell(i_Cell,:)                       = CI_speed;
            CI_speed_90_cell(i_Cell,:)                    = CI_speed_90;
            R_boot_mean_speed_cell(i_Cell,1)              = R_boot_mean_speed;
            HalvesCorr_Pearson_Boot_cell_shuf{i_Cell,1}     = HalvesCorr_Pearson_shuf;
            p_rank_ratemap_cell(i_Cell,1)                   = p_rank;
            Cliff_delta_ratemap_cell(i_Cell,1)              = Cliff_delta;
            HalvesCorr_speed_Pearson_Boot_cell_speed{i_Cell,1}  = HalvesCorr_speed_Pearson_shuf;
            p_rank_speed_cell(i_Cell,1)                         = p_rank_speed;
            Cliff_delta_speed_cell(i_Cell,1)                    = Cliff_delta_speed;


            % set output folder

            OutDir_Boot = fullfile(saveFolder, 'RateMapBootstrap_shuffle_260204');
            [status,msg,msgID] = mkdir(OutDir_Boot);

            animalID = sprintf('Mouse%sDay%d', mice_str_Folder{mice_ind}, Day);
            cellID   = sprintf('CellID#%d', i_Cell);
            fname = sprintf('RateMapBootstrap %s_%s.jpeg', animalID, cellID);
            exportgraphics(fig_boot, fullfile(OutDir_Boot, fname));


            overallTime_perCell = toc(overallTic);
            fprintf('Total time per cell: %.3f s\n', overallTime_perCell);

        end

        %% summarize some outputs into table
        %% save variables

        % SpkDetectionResults.CaProcessingParm = CaProcessingParm; % saving raw traces is heavy
        SpkDetectionResults_save.SpkTime = SpkTime;
        SpkDetectionResults_save.Peaks = Peaks;
        SpkDetectionResults_save.SpkTimeAligned_table = SpkTimeAligned_table;
        SpkDetectionResults_save.FiringRate = FiringRate_cell;

        RateMaps.SpkDetection.runSpkIdx = runSpkIdx_sp;
        RateMaps.SpkDetection.neuronPos = neuronPos_sp;
        RateMaps.SpkDetection.numEvent = numEvent_sp;
        RateMaps.SpkDetection.SpeedResults = SpeedResults;

        RateMaps.Adapt.speedInd_list = speedInd_list;
        RateMaps.Adapt.binSize_cm = binSize_cm;
        RateMaps.Adapt.rateMap_ori    = rateMap_ori_sp;
        RateMaps.Adapt.rateMap_thr    = rateMap_thr_sp;
        RateMaps.Adapt.occMap_ori     = occMap_ori_sp;
        RateMaps.Adapt.occMap_thr     = occMap_thr_sp;
        RateMaps.Adapt.spikeMap_Rawcount = spikeMap_Rawcount_sp;
        RateMaps.Adapt.areaMap = areaMap_sp;
        RateMaps.Adapt.XYEdges = Edges_sp;
        RateMaps.Adapt.Dur_Running = Dur_Running_sp;

        RateMaps.Uniform.XY_oldnew_cell             = XY_oldnew_cell;
        RateMaps.Uniform.occMapNormalized           = occMapNormalized_sp;
        RateMaps.Uniform.occMapUniform_Gaussian     = occMapUniform_Gaussian_sp;
        RateMaps.Uniform.rateMapUniform             = rateMapUniform_sp;
        RateMaps.Uniform.rateMapUniform_Gaussian    = rateMapUniform_Gaussian_sp;

        RateMaps.FineRateMap.rateMap_Fine_sp = rateMap_Fine_sp;
        RateMaps.FineRateMap.PeakFiringRate_sp = PeakFiringRate_sp;
        RateMaps.FineRateMap.xyFine_sp = [xs_sp, ys_sp];
        RateMaps.FineRateMap.factor = factor;


        TotalRecordingTime = Time_aligned(end);

        if count == n_Cell % to avoid overwriting
            %% save basics
            s_savemat = strcat('SpaceAnalysis_', mice_str{mice_ind}, '_d', num2str(Day));
            ss = fullfile(saveFolder_Basic, s_savemat);
            save(ss, 'SpkDetectionResults_save', 'RateMaps',  'FieldAnalysis_cell', 'SpInfo_cell', 'Coherence_cell', ...
                'SpInfo_shuf_cell', 'Coherence_shuf_cell', 'speedModResult_cell', ...
                'HalveCompareField_stability_cell', 'HalveCompareSpeed_stability_cell', 'HalveCompareField_reliability_cell', 'HalveCompareSpeed_reliability_cell',...
                'Halves_shuf_stab_cell', 'Halves_shuf_reliab_cell', ...
                'VideoSampling', 'TotalRecordingTime', ...
                '-v7.3');

            %% save field analysis results
            s_savemat = strcat('FieldAnalysisRedo_', mice_str{mice_ind}, '_d', num2str(Day));
            ss = fullfile(OutDir_FRD, s_savemat);
            save(ss, 'FRD_FieldAnalysis_Redo_cell', 'FRD_ACGanalysis_cell', ...
                '-v7.3');

            %% save bootstrap results
            s_savemat = strcat('RateMapBootstrap_Shuffle', mice_str{mice_ind}, '_d', num2str(Day));
            ss = fullfile(OutDir_Boot, s_savemat);

            save(ss, 'HalvesCorr_Pearson_even_cell', 'HalvesCorr_speed_Pearson_even_cell', ...
                "HalvesCorr_Pearson_Boot_cell", "CI_rateMap_cell", "R_boot_mean_rateMap_cell",...
                "HalvesCorr_speed_Pearson_Boot_cell", "CI_speed_cell", "R_boot_mean_speed_cell",...
                "CI_rateMap_90_cell", 'CI_speed_90_cell', ...
                'HalvesCorr_Pearson_Boot_cell_shuf', 'p_rank_ratemap_cell', 'Cliff_delta_ratemap_cell', ...
                'HalvesCorr_speed_Pearson_Boot_cell_speed', 'p_rank_speed_cell', 'Cliff_delta_speed_cell',...
                '-v7.3');

        end


    end
end




%% functions
%% Spike Detection
function [CaProcessingParm, SpkData, Peaks] = func_SpkDetection(TraceData, thr)
%%

% Step 0: Ensure ImgData is numeric
if istable(TraceData)
    TraceData = table2array(TraceData);
end
if ~isnumeric(TraceData)
    error('Input ImgData must be numeric or convertible to numeric.');
end

% Step 1: Extract time vector
SpikeTimes = TraceData(:,1);
TraceData(:,1) = [];

% Step 2: Clean NaNs
rowsToRemove = find(isnan(TraceData(:,1)));
if ~isempty(rowsToRemove)
    TraceData(rowsToRemove, :) = [];
    SpikeTimes(rowsToRemove, :) = [];
end
cleanedTraceData = TraceData;

% Step 3: Initialize
[numRows, numTraces] = size(TraceData); %#ok<ASGLU>
threshold = zeros(1, numTraces);
baseline  = zeros(1, numTraces);
prominence_used = zeros(1, numTraces);
SpkData   = cell(1, numTraces);
Peaks = cell(1, numTraces);

% std_dF = std(TraceData, 0, 1);
MAD = mad(TraceData, 1, 1);

% Step 4: For each ROI
for idata = 1:numTraces
    currentData = TraceData(:, idata);

    % Baseline and threshold
    baseline(idata) = 0; % not used, already corrected
    threshold(idata) = (thr * MAD(idata));

    % Prominence (can tune this factor if needed)
    prominence_used(idata) = 0; % not used


    Fs = size(TraceData,1)/SpikeTimes(end,1);
    [pks, locs2] = findpeaks(currentData, 'MinPeakHeight', threshold(idata), ...
        'MinPeakProminence', threshold(idata)/2, 'MinPeakDistance', 0.1 * Fs);
    %'MinPeakProminence', prominence_used(idata));

    SpkData{idata} = SpikeTimes(locs2);
    Peaks{idata} = pks;
end

%% Output
CaProcessingParm = struct(...
    'threshold', threshold, ...
    'baseline', baseline, ...
    'prominence', prominence_used, ...
    'cleanedTraceData', cleanedTraceData, ...
    'SpikeTimes', SpikeTimes ...
    );
end


%% Align spikes to video timestamps
function [SpkTimeAligned, SpkTimeAligned_table] = func_AlignSpikeToVideoTime(Table_video, SpkTime)
%%
timeAligned = Table_video.time_aligned(:); %read behavior video timestamps which was aligned to Ca imaging recording start
dsArray = table2array(Table_video);
varNames    = Table_video.Properties.VariableNames;
nCols  = size(dsArray, 2);

numCells = numel(SpkTime);
SpkTimeAligned = cell(1, numCells);
SpkTimeAligned_table = cell(1, numCells);

for iCell = 1:numCells
    spikeTimes = SpkTime{iCell}(:);
    nSpikes = numel(spikeTimes);

    % Find nearest induces using interp1
    % this calculation does not lose any spike data because Fs of animal behavior
    % video is faster than Ca imaging
    idx = interp1(timeAligned, 1:numel(timeAligned), spikeTimes, 'nearest', NaN);

    % If the nearest time difference is more than 1 second, assign NaN row
    mappedData = nan(nSpikes, nCols);
    valid = ~isnan(idx);
    if any(valid)
        diffT = abs(timeAligned(idx(valid)) - spikeTimes(valid));
        keep   = diffT <= 1;
        rows   = idx(valid(keep));
        mappedData(valid(keep), :) = dsArray(rows, :);
    end

    SpkTimeAligned{iCell} = mappedData;
    SpkTimeAligned_table{iCell} = array2table(mappedData, 'VariableNames', varNames);
end


end



%% computeSpatialRateMap

% gaussian filter ignoring NaN
function Zs = func_gaussSmoothIgnoreNaN(Z, sigma, kernelSize) %#ok<DEFNU>
% Z … occupancy or spike map (with NaN)
% sigma … gaussian σ (pixels)

% create gaussian kernel
hsize = kernelSize(1);

if mod(hsize,2)==0, hsize=hsize+1; end
g = fspecial('gaussian', hsize, sigma);

W = ~isnan(Z); %prepare mask and data
Z(~W) = 0;

num = imfilter(Z, g, 'replicate', 'conv'); %convolve data and mask
den = imfilter(double(W), g, 'replicate', 'conv');

Zs = num ./ den; % normalizing
Zs(den==0) = NaN;
end


%% Adaptive bin rate map

function [occMapMerged, spkMapMerged, rateMapMerged, xEdges, yEdges, nBin, Dur_Running, areaMap, rateMapMerged_ori, occMapMerged_ori] = ...
    func_adaptiveBinRateMap(mousePos_run, neuronPos, AnimalTrack, binSize_cm, VideoSampling, occThresh, sigma, kernelSize) %#ok<INUSD>
%% initialize

minOccSec = occThresh;
pos = mousePos_run;

cropRect = AnimalTrack.TrackData.cropRect;
psh = AnimalTrack.TrackData.psh;
psw = AnimalTrack.TrackData.psw;
xDim = cropRect(3)*psw;
yDim = cropRect(4)*psh;
AreanaxEdges = 0:binSize_cm:ceil(xDim/ binSize_cm)*binSize_cm;
AreanayEdges = 0:binSize_cm:ceil(yDim / binSize_cm)*binSize_cm;

A = max(pos(:,1)) - min(pos(:,1));
B = max(pos(:,2)) - min(pos(:,2));

N = round(max([A,B]) / binSize_cm);
nBin = N;

% ==== 1. Equal-occupancy binning ====
if ~isempty(pos)
    nPos = size(pos,1);
    xSorted = sort(pos(:,1));
    ySorted = sort(pos(:,2));
    idxX = round(linspace(1, nPos, nBin+1));
    idxY = round(linspace(1, nPos, nBin+1));
    xEdges = unique(xSorted(idxX));
    yEdges = unique(ySorted(idxY));
else
    xEdges = (0:binSize_cm:AreanaxEdges(end))';
    yEdges = (0:binSize_cm:AreanayEdges(end))';
end

xEdges = [0; xEdges; AreanaxEdges(end); AreanaxEdges(end)+1];
yEdges = [0; yEdges; AreanayEdges(end); AreanayEdges(end)+1];
xEdges = unique(xEdges);
yEdges = unique(yEdges);


% ==== 2. initialize map ====
occCount = histcounts2(pos(:,1), pos(:,2), xEdges, yEdges)';
spkCount = histcounts2(neuronPos(:,1), neuronPos(:,2), xEdges, yEdges)';
occMap = occCount / VideoSampling; % [count] to [time]


% ==== 3. merge bins with < 2sec threshold ====
occMapMerged = occMap;
spkMapMerged = spkCount;
[nRows, nCols] = size(occMapMerged);

% merge from center
maskLow = occMapMerged < minOccSec;
while any(maskLow(:))
    [rLow, cLow] = find(maskLow);

    % priotize center bins
    [~, idxSort] = sortrows(abs([rLow cLow] - [nRows nCols]/2));
    rLow = rLow(idxSort); cLow = cLow(idxSort);

    for k = 1:length(rLow)
        r = rLow(k); c = cLow(k);
        if isnan(occMapMerged(r,c)) || occMapMerged(r,c) >= minOccSec
            continue; % finished marging
        end

        % neighboring bin coordinates
        neighR = max(r-1,1):min(r+1,nRows);
        neighC = max(c-1,1):min(c+1,nCols);

        % remove self bin
        neighOcc = occMapMerged(neighR, neighC);
        neighOcc(r - neighR(1)+1, c - neighC(1)+1) = -Inf;

        % chose max occupancy bin
        [~, idxMax] = max(neighOcc(:));
        [nr, nc] = ind2sub(size(neighOcc), idxMax);
        targetR = neighR(nr);
        targetC = neighC(nc);

        % merge
        occMapMerged(targetR, targetC) = occMapMerged(targetR, targetC) + occMapMerged(r,c);
        spkMapMerged(targetR, targetC) = spkMapMerged(targetR, targetC) + spkMapMerged(r,c);

        % change original bin to NaN
        occMapMerged(r,c) = NaN;
        spkMapMerged(r,c) = NaN;
    end
    maskLow = occMapMerged < minOccSec;
end
rateMapMerged = spkMapMerged ./ occMapMerged;


% are map
[nRows, nCols] = size(occMapMerged);
binWidth  = diff(xEdges);
binHeight = diff(yEdges);
areaMap = nan(nRows, nCols);
for r = 1:nRows
    for c = 1:nCols
        if isnan(occMapMerged(r,c))
            continue
        end
        areaMap(r,c) = binWidth(c) * binHeight(r); % cm^2
    end
end

% normalize occupancy
occMapNormalized = occMapMerged ./ areaMap;  % sec/cm^2


% thresholding by occupancy density and bin area size
rateMapMerged_ori = rateMapMerged;
occMapMerged_ori = occMapMerged;
Dur_Running = sum(occMap(:), 'omitmissing');

thr_area = 50;
rateMapMerged(areaMap > thr_area) = 0;
occMapMerged(areaMap > thr_area) = 0;

thr_sec_perarea = 0.1;
rateMapMerged(occMapNormalized < thr_sec_perarea) = 0;
occMapMerged(occMapNormalized < thr_sec_perarea) = 0;

end



%% Generate and Show Gaussian Occ and RateMap
function [rateMapUniform, rateMapUniform_Gaussian, occMapNormalized, occMapGaussian, X_old, Y_old, X_new, Y_new, occMapUniform_Gaussian] = ...
    func_MakeAndShowGaussianRateMap(AnimalTrack, occMap_ori_sp, Edges_sp, thr_sp, i_Cell, areaMap_sp, binSize_cm, thresholds_speed, sigma, kernelSize, rateMap_thr_sp, doShow)

%%
occMapMerged = occMap_ori_sp{thr_sp,1};
Edges = Edges_sp{thr_sp,i_Cell};
areaMap = areaMap_sp{thr_sp,i_Cell};
rateMapMerged = rateMap_thr_sp{thr_sp,i_Cell};

[occMapNormalized, occMapGaussian, occMapUniform, occMapUniform_Gaussian, X_old, Y_old, X_new, Y_new, xPlot, yPlot, xDim, yDim] = ...
    computeGaussianOccMap(AnimalTrack, occMapMerged, Edges, areaMap, binSize_cm, sigma, kernelSize);
[rateMapUniform, rateMapUniform_Gaussian] = computeGaussianRateMap(rateMapMerged, X_old, Y_old, X_new, Y_new, sigma, kernelSize);

% plot maps
if doShow
    plotGaussianRateMaps(occMapNormalized, occMapGaussian, rateMapUniform_Gaussian, ...
        xPlot, yPlot, X_new, Y_new, xDim, yDim, thr_sp, thresholds_speed, occMapUniform, occMapUniform_Gaussian);
end

end

% ---- helper functions ----
% ----- compute uniform and gaussian occ map
function [occMapNormalized, occMapGaussian, occMapUniform, occMapUniform_Gaussian, X_old, Y_old, X_new, Y_new, xPlot, yPlot, xDim, yDim] = ...
    computeGaussianOccMap(AnimalTrack, occMapMerged, Edges, areaMap, binSize_cm, sigma, kernelSize)

cropRect = AnimalTrack.TrackData.cropRect;
psh = AnimalTrack.TrackData.psh; psw = AnimalTrack.TrackData.psw;
xDim = cropRect(3)*psw; yDim = cropRect(4)*psh;

xEdges = Edges{1};
yEdges = Edges{2};
xPlot = xEdges(1:end-1);
yPlot = yEdges(1:end-1);

% normalize occ map
occMapNormalized = occMapMerged ./ areaMap;
occMapNormalized(isnan(occMapNormalized)) = 0;

% make grid
[X_old, Y_old] = meshgrid((xEdges(1:end-1)+xEdges(2:end))/2, ...
    (yEdges(1:end-1)+yEdges(2:end))/2);
xGrid = 0:binSize_cm:ceil(xDim/binSize_cm)*binSize_cm;
yGrid = 0:binSize_cm:ceil(yDim/binSize_cm)*binSize_cm;
[X_new, Y_new] = meshgrid(xGrid, yGrid);

% smoothing occ map
occMapGaussian = imgaussfilt(occMapNormalized, sigma, ...
    "FilterSize", kernelSize, 'Padding', 'symmetric');

mask = isnan(occMapMerged);
occMapMerged(mask) = 0;
occMapUniform = interp2(X_old, Y_old, occMapMerged, X_new, Y_new, 'linear'); % seconds
occMapUniform_Gaussian = interp2(X_old, Y_old, occMapGaussian, X_new, Y_new, 'linear'); % seconds

end

% ----- compute uniform and gaussian rate map
function [rateMapUniform, rateMapUniform_Gaussian] = computeGaussianRateMap(rateMapMerged, X_old, Y_old, X_new, Y_new, sigma, kernelSize)
rateMapMerged(isnan(rateMapMerged)) = 0;

rateMapUniform = interp2(X_old, Y_old, rateMapMerged, X_new, Y_new, 'linear'); % interpolate to uniform grid
rateMapUniform(isnan(rateMapUniform)) = 0;

rateMapUniform_Gaussian = imgaussfilt(rateMapUniform, sigma, 'FilterSize', kernelSize, 'Padding', 'replicate'); % Gaussian smoothing
end

% ----- plot rate maps
function plotGaussianRateMaps(occMapNormalized, occMapGaussian, rateMap_Gaussian, xPlot, yPlot, X_new, Y_new, xDim, yDim, thr_sp, thresholds_speed, occMapUniform, occMapUniform_Gaussian) %#ok<INUSD>
%%
colormap jet

subplot(5,6,12 + thr_sp)
pcolor(xPlot, yPlot, occMapNormalized); shading flat;
axis tight equal off; set(gca, 'YDir','reverse')
xlim([0 xDim]); ylim([0 yDim]);
title(sprintf('Occupancy map at %d cm/s', thresholds_speed(thr_sp)))

subplot(5,6,18 + thr_sp)
pcolor(X_new, Y_new, occMapUniform_Gaussian); shading flat;
axis tight equal off; set(gca, 'YDir','reverse')
xlim([0 xDim]); ylim([0 yDim]);
title('Gaussian filtered occ map')

subplot(5,6,24 + thr_sp)
pcolor(X_new, Y_new, rateMap_Gaussian); shading flat;
axis tight equal off; set(gca, 'YDir','reverse')
xlim([0 xDim]); ylim([0 yDim]);
title('Gaussian filtered rate map')

end


%% show CaTrace and Spikes
function [fig, FiringRate] = func_showCaTraceAndSpikes_251013(AnimalTrack, Table_video, SpkTimeAligned_table, TraceData, runSpkIdx_sp, i_Cell, CaProcessingParm, SpeedResults)
% arena size
cropRect = AnimalTrack.TrackData.cropRect;
psh = AnimalTrack.TrackData.psh; psw = AnimalTrack.TrackData.psw;
xDim = cropRect(3)*psw; yDim = cropRect(4)*psh;
xmin = 0; xmax = xDim; ymin = 0; ymax = yDim;

% spike positions
x_all = Table_video.x; y_all = Table_video.y;
spikeTime = SpkTimeAligned_table{1,i_Cell};
x_spk = spikeTime.x; y_spk = spikeTime.y;

% figure
speedInd_list = [1 2 3 4 5 7 11 16 21];
fig = figure('Position',[-1253 73 1120 876],'Visible','off');
% fig = figure('Position',[-1253 73 1120 876],'Visible','on');

% --- Calcium trace ---
subplot(5,6,[1 2])
plot(TraceData(:,1), TraceData(:,i_Cell+1),'k'), axis tight, box off
Fs = size(TraceData,1)/TraceData(end,1);
[NonRecordingDur, ~] = func_detectNonRecording(TraceData(:,2), Fs);
FiringRate = length(x_spk)/(TraceData(end,1) - NonRecordingDur);
title(sprintf('Firing rate = %.3f Hz; %d events', FiringRate, numel(x_spk)))
hold on, yline(CaProcessingParm.threshold(i_Cell),'--','Color',[1 0 0 0.6])
xlabel('time [s]')
ylabel('Calcium ΔF/F')

v = SpeedResults.smoothedSpeed;
t = Table_video.time_aligned;
y_range = max(TraceData(:,i_Cell+1));
v_scaled = v / (max(v)) * 0.2 * y_range;  %
v_scaled = v_scaled + 1.1*y_range;  % plot above trace
plot(t, v_scaled, 'Color', [0.3, 0.7, 1, 0.6])

ax = gca;
yyaxis left
if max(v_scaled)>0; ylim([0 max(v_scaled)]); end
ax.YColor = 'k';

yyaxis right
if max(v_scaled)>0; ylim([0 max(v_scaled)]); end
yticks([1.1*y_range, max(v_scaled)])
yticklabels({ '0', sprintf('%.1f', max(v)) })
ylabel('Speed [cm/s]')

% ax.YAxis(2).Color = 'none';
ax.YColor = [0.3, 0.7, 1];
ax.YLabel.Color = [0.3, 0.7, 1];

% --- rate maps with speed thresholds ---
thrVals = SpeedResults.Thresholds(speedInd_list);
for k = 1:5
    spId   = speedInd_list(k);
    runIdx = runSpkIdx_sp{spId,i_Cell};

    % spike scatter
    subplot(5,6,6+k); hold on;
    plotScatter(x_all,y_all,x_spk(runIdx),y_spk(runIdx),xmin,xmax,ymin,ymax,thrVals(k),nnz(runIdx))
end
end

% --- helper functions ---
function plotScatter(x_all,y_all,x_spk,y_spk,xmin,xmax,ymin,ymax,thrVal,numEvent)
plot(x_all,y_all,'Color',[0.8 0.8 0.8],'LineWidth',1), hold on
scatter(x_spk,y_spk,10,'r','filled')
if xmin == xmax; xmax = xmin+1; end
if ymin == ymax; ymax = ymin+1; end
axis equal off, xlim([xmin xmax]), ylim([ymin ymax]), set(gca,'YDir','reverse')
title(sprintf('Spike at > %.1f cm/s, %d events', thrVal,numEvent))
end

function [NonRecordingDur, idx_remove] = func_detectNonRecording(trace, Fs)
tol = 1e-12;
isConst = [false; abs(diff(trace)) <= tol];
grp = bwlabel(isConst);
stats = regionprops(grp,'Area'); %#ok<MRPBW>
minLen = 100;
idx_remove = ismember(grp, find([stats.Area] >= minLen));
NonRecordingDur = nnz(idx_remove)/Fs;
end

%% Place field size and other parameters
function [fig_field, FieldStats_ratio, FieldStats_STD] = func_FieldAnalysis(rateMap_Fine_sp, thr_sp, binSize_cm, factor)
%%
rateSmooth = rateMap_Fine_sp{thr_sp};
Length_per_bin = binSize_cm/factor; % cm
Area_per_bin   = Length_per_bin^2;

rateSmooth_NaN = rateSmooth;
rateSmooth_NaN(rateSmooth_NaN == 0) = NaN;

fig_field = [];

% --- ratio thresholding ---
thr_field_ratio_list = 0.0:0.2:0.8;
FieldStats_ratio = func_runFieldAnalysis(rateSmooth, rateSmooth_NaN, ...
    thr_field_ratio_list, Area_per_bin, "ratio");

% --- MAD thresholding ---
% thr_field_STD_list = [2 3 4 5 6];
thr_field_STD_list = [0 1 1.5 2 3];
FieldStats_STD = func_runFieldAnalysis(rateSmooth, rateSmooth_NaN, ...
    thr_field_STD_list, Area_per_bin, "mad");
end

% --- helper functions ---
function FieldStats = func_runFieldAnalysis(rateSmooth, rateSmooth_NaN, thr_list, Area_per_bin, mode)
nThr = numel(thr_list);
FieldStats = cell(nThr, 12);

for iThr = 1:nThr
    % --- threshold
    if mode == "ratio"
        peakRate = max(rateSmooth(:), [], 'omitnan');
        thr_val  = thr_list(iThr) * peakRate;
    elseif mode == "mad" % MAD
        thr_val = std(rateSmooth(:)) * thr_list(iThr); % STDEV(σ) ≈ 1.25*MAD
    end

    [rateMap_thr, stats, inFR, outFR] = func_thresholdAndExtract(rateSmooth, rateSmooth_NaN, thr_val, Area_per_bin);
    [~, fieldScores, allFieldsScore, Dist] = func_ComputeBorderScore(rateMap_thr, thr_val);
    [PeakFR, FieldSize_Peak, FieldSize_All, FieldNum] = func_summarizeStats_Field(stats);
    FieldStats(iThr,:) = {thr_list(iThr), rateMap_thr, stats, inFR, outFR, ...
        PeakFR, FieldSize_Peak, FieldSize_All, FieldNum, fieldScores, Dist, allFieldsScore};

end

FieldStats = cell2table(FieldStats, ...
    "VariableNames",["FDthreshold","thresholded_image","stats","infield_meanfiring","outfield_meanfiring",...
    "PeakFR","FieldSize_Peak","FieldSize_All","FieldNum", "fieldScores", "DistanceFromBorder_px", "allFieldsScore"]);
end

function [PeakFR, FieldSize_Peak, FieldSize_All, FieldNum] = func_summarizeStats_Field(stats)
if ~isempty(stats)
    [PeakFR,I]   = max([stats.MaxIntensity]);
    FieldSize_Peak = stats(I).Area_real;
    FieldNum     = numel(stats);
    FieldSize_All = sum([stats.Area_real]);
else
    PeakFR = 0; FieldSize_Peak = 0; FieldNum = 0; FieldSize_All = 0;
end
end

function func_drawFieldMap(rateMap_thr, iThr, mode, thr, PeakFR, FieldSize_Peak, FieldSize_All, FieldNum, borderScore, allFieldsScore) %#ok<DEFNU>
if mode == "ratio"
    nexttile(iThr);
    imagesc(rateMap_thr); axis tight equal off; colormap jet
    title(sprintf('thr: %.1f * peak', thr));
else
    nexttile(5+iThr);
    imagesc(rateMap_thr); axis tight equal off; colormap jet
    title(sprintf('thr: %.1f * STD', thr));

    nexttile(10+iThr);
    s = {sprintf('PeakFR = %.3f Hz', PeakFR), ...
        sprintf('FieldSizePeak = %.0f cm2', FieldSize_Peak), ...
        sprintf('FieldSizeAll = %.0f cm2', FieldSize_All), ...
        sprintf('FieldNum = %d', FieldNum), ...
        sprintf('BorderScore = %.2f', borderScore), ...
        sprintf('AllFieldsScore = %.2f', allFieldsScore)};
    text(0,1,s, 'VerticalAlignment','top','Units','normalized');
    axis off
end
end

function [rateMap_thr, stats, inFR, outFR] = func_thresholdAndExtract(rateMap, rateMapNaN, thr, Area_per_bin)
mask = rateMap <= thr;
rateMap_thr = rateMap;
rateMap_thr(mask) = 0;

BW = ~mask;
stats = regionprops(BW, rateMap_thr, ...
    "Area","Centroid","Circularity","MajorAxisLength","MinorAxisLength", ...
    "Orientation","MaxIntensity","MeanIntensity","WeightedCentroid","PixelIdxList");

for i = 1:numel(stats)
    [~, idxLocal] = max(rateMap_thr(stats(i).PixelIdxList));
    maxIdx = stats(i).PixelIdxList(idxLocal);
    [row,col] = ind2sub(size(rateMap_thr), maxIdx);
    stats(i).MaxIPos = [col,row];
    stats(i).Area_real = stats(i).Area * Area_per_bin;
end

% in/out field firing rates
inFieldValues  = rateMapNaN(~mask);
outFieldValues = rateMapNaN(mask);
inFR  = mean(inFieldValues,'omitnan');
outFR = mean(outFieldValues,'omitnan');
end


%% border score
function [borderScore, fieldScores, allFieldsScore, Dist] = func_ComputeBorderScore(rateMap, thr)
% https://www.science.org/doi/10.1126/science.1166466#supplementary-materials
% https://www.nature.com/articles/s41593-017-0055-3#Sec11
% ‘Border cells’ were defined as cells with border scores above 0.5.

% ComputeBorderScore - calculate border score from a 2D rate map
%
% INPUTS:
%   rateMap : 2D matrix (NaN for unvisited bins, numeric for firing rates)
%   thr     : threshold for field detection (default = 0.2 * peak)
%
% OUTPUTS:
%   borderScore : max border score among detected fields
%   fieldScores : border score for each individual field

if nargin < 2
    thr = 0.2 * max(rateMap(:)); % 20% peak threshold
end

% --- thresholding to detect fields ---
fieldMask = rateMap >= thr;
CC = bwconncomp(fieldMask); % connected components
stats = regionprops(CC, 'PixelIdxList', 'Centroid');

if isempty(stats)
    borderScore = NaN;
    fieldScores = [];
    allFieldsScore = NaN;
    Dist = NaN;
    return;
end

[nRows, nCols] = size(rateMap);
fieldScores = nan(1, numel(stats));
Dist = nan(1, numel(stats));
for k = 1:numel(stats)
    idx = stats(k).PixelIdxList;
    [r, c] = ind2sub([nRows, nCols], idx);

    % --- centroid-to-wall distance (dc) ---
    centroid = stats(k).Centroid;
    dc = min([centroid(1)-1, nCols-centroid(1), centroid(2)-1, nRows-centroid(2)]);

    % --- fraction of field pixels touching wall (df) ---
    onWall = (r==1 | r==nRows | c==1 | c==nCols);
    df = sum(onWall) / numel(idx);

    % --- border score ---
    fieldScores(k) = (df - dc/(max(nRows,nCols))) / (df + dc/(max(nRows,nCols)));

    Dist(k) = dc; %distance from border to center
end

% neuron border score = max over fields
borderScore = max(fieldScores);

% --- all-fields combined score ---
allMask = fieldMask;  % all fields merged
[r, c] = ind2sub([nRows, nCols], find(allMask));
centroidAll = mean([c, r], 1);
dcAll = min([centroidAll(1)-1, nCols-centroidAll(1), centroidAll(2)-1, nRows-centroidAll(2)]) / max(nRows, nCols);

onWallAll = (r==1 | r==nRows | c==1 | c==nCols);
dfAll = sum(onWallAll) / numel(r);

allFieldsScore = (dfAll - dcAll) / (dfAll + dcAll);

end

%% Spatial Info
function [meanRate, info_bits_per_spike, info_bits_per_sec, spatialSparsity, spatialSelectivity] = func_SkaggsSpatialInfo_251014(rateMap, occMap)
% SkaggsSpatialInfo  compute Skaggs spatial information (bits/spike)
% Inputs:
%   rateMap : MxN matrix of firing rates [Hz] for each spatial bin
%   occMap  : MxN matrix of occupancy time [s] in each spatial bin
%
% Outputs:
%   meanRate : occupancy-weighted mean firing rate [Hz]
%   info_bits_per_spike : Skaggs information in bits/spike
%   info_bits_per_sec   : bits/sec = info_bits_per_spike * meanRate
%%

r   = rateMap(:);
occ = occMap(:);

valid = ~isnan(r) & ~isnan(occ); % & (occ > occThresh); % only bins with enough occupancy
if ~any(valid)
    meanRate = NaN;
    info_bits_per_spike = NaN;
    info_bits_per_sec = NaN;
    spatialSparsity = NaN;
    return;
end

% total time and mean rate (occupancy-weighted)
T = sum(occ(valid));             % total occupied time (s)
meanRate = sum(r(valid) .* occ(valid)) / T;  % Hz

% avoid degenerate case
if meanRate <= 0
    info_bits_per_spike = NaN;
    info_bits_per_sec = NaN;
    spatialSparsity = NaN;
    spatialSelectivity = NaN;
    return;
end

% Skaggs spatial sparsity
meanRateSq = sum(occ(valid) .* r(valid).^2) / T;
spatialSparsity = (meanRate^2) / meanRateSq;

spatialSelectivity = 1 - spatialSparsity;

% Skaggs spatial information (bits/spike)
% compute p_i and lambda_i/meanRate only for bins with lambda>0
lambda = r(valid);
p_i = occ(valid) / T;            % occupancy probability

positive = lambda > 0;           % include only lambda>0 (0 term -> 0 in limit)
ratio = lambda(positive) / meanRate;  % lambda_i / lambda_bar
p_pos = p_i(positive);

info_bits_per_spike = sum( p_pos .* ratio .* log2(ratio) ); % Skaggs formula (bits per spike)
info_bits_per_sec = info_bits_per_spike * meanRate; % bits per second
end

%% shuffling for spatial info
function [SpInfo_shuf, Coherence_shuf, shuffIdx, neuronPos_shuff] = ...
    func_ShuffleTest_adapt(thr_sp, SpeedResults, spikeFrames, mousePos, AnimalTrack, occThresh, sigma, kernelSize, VideoSampling, binSize_cm, SpInfo, Coherence, spikePos)
%%
info_perspike   = SpInfo.info_perspike;
info_persec     = SpInfo.info_persec;
spatialSparsity = SpInfo.spatialSparsity;
spatialSelectivity = SpInfo.spatialSelectivity;
Z = Coherence.Z;

%% spike index during running
nShuffle = 1000;
runFrames = SpeedResults.idx{thr_sp};

[Idx_Lia, ~] = ismember(runFrames, spikeFrames); %spike index during running
nRunFrames = length(Idx_Lia);


%% generate occ map
runFrames = SpeedResults.idx{thr_sp};
runSpkIdx = ismember(spikeFrames, runFrames);
neuronPos = spikePos(runSpkIdx, :);
mousePos_run = mousePos(runFrames,:);

[occMapMerged, ~, ~, xEdges, yEdges, ~, ~, areaMap, ~, ~] = ...
    func_adaptiveBinRateMap(mousePos_run, neuronPos, AnimalTrack, binSize_cm, VideoSampling, occThresh, sigma, kernelSize);
Edges = {xEdges, yEdges};
[~, ~, ~, occMapUniform_Gaussian, X_old, Y_old, X_new, Y_new, ~, ~, ~, ~] = ...
    computeGaussianOccMap(AnimalTrack, occMapMerged, Edges, areaMap, binSize_cm, sigma, kernelSize);

% figure
% subplot(121)
% imagesc(occMapUniform_Gaussian)
% subplot(122)
% imagesc(rateMapUniform_Gaussian)
% colormap jet

%% shuffling

info_shuf = nan(nShuffle,1);
meanRate_shuf = nan(nShuffle,1);
info_sec_shuf = nan(nShuffle,1);
spatialSparsity_shuf = nan(nShuffle,1);

shuffIdx = cell(nShuffle,1);
neuronPos_shuff= cell(nShuffle,1);

Z_shuf = nan(nShuffle,1);
Coherence_SpCorr = nan(nShuffle,1);

% tic
for i_shuffle = 1:nShuffle %usually for is faster than parfor
    % parfor i_shuffle = 1:nShuffle
    shift = randi(nRunFrames);
    local_posRun = mousePos_run;

    shuffIdx{i_shuffle} = circshift(Idx_Lia, shift); % spike timing is shifted within running duration
    neuronPos_shuff{i_shuffle} = local_posRun(shuffIdx{i_shuffle}, :);


    [~, ~, rMap_shuf, ~, ~, ~, ~, ~, ~, ~] = ...
        func_adaptiveBinRateMap(mousePos_run, neuronPos_shuff{i_shuffle}, AnimalTrack, binSize_cm, VideoSampling, occThresh, sigma, kernelSize);
    [rateMapUniform_shuf, rateMapUniform_Gaussian_shuf] = computeGaussianRateMap(rMap_shuf, X_old, Y_old, X_new, Y_new, sigma, kernelSize);

    [meanRate_shuf(i_shuffle), info_shuf(i_shuffle), info_sec_shuf(i_shuffle), spatialSparsity_shuf(i_shuffle)] = func_SkaggsSpatialInfo_251014(rateMapUniform_Gaussian_shuf, occMapUniform_Gaussian);

    % [Z_shuf(i_shuffle), Coherence_SpCorr(i_shuffle), ~] = func_SpatialCoherence(rateMapUniform_Gaussian_shuf);
    [Z_shuf(i_shuffle), Coherence_SpCorr(i_shuffle), ~] = func_SpatialCoherence_251014(rateMapUniform_shuf);
end
% toc %2sec

spatialSelectivity_shuf = 1 - spatialSparsity_shuf;


% p-test (one-side)
info_obs = info_perspike;
pval_spinfo     = (sum(info_shuf >= info_obs) + 1) / (nShuffle + 1);
pval_spsecinfo  = (sum(info_sec_shuf >= info_persec) + 1) / (nShuffle + 1);
pval_sparsity   = (sum(spatialSparsity_shuf <= spatialSparsity) + 1) / (nShuffle + 1);
pval_selectivity   = (sum(spatialSelectivity_shuf >= spatialSelectivity) + 1) / (nShuffle + 1);
pval_coherence  = (sum(Z_shuf >= Z) + 1) / (nShuffle + 1);

% return outputs
SpInfo_shuf.info_obs        = info_obs;
SpInfo_shuf.info_shuf       = info_shuf;
SpInfo_shuf.meanRate_shuf   = meanRate_shuf;
SpInfo_shuf.info_obs_sec        = info_persec;
SpInfo_shuf.info_sec_shuf   = info_sec_shuf;
SpInfo_shuf.pval_spinfo     = pval_spinfo;
SpInfo_shuf.pval_spsecinfo  = pval_spsecinfo;
SpInfo_shuf.pval_sparsity   = pval_sparsity;
SpInfo_shuf.spatialSparsity_shuf   = spatialSparsity_shuf;
SpInfo_shuf.pval_selectivity   = pval_selectivity;

Coherence_shuf.Z_shuf   = Z_shuf;
Coherence_shuf.Coherence_SpCorr   = Coherence_SpCorr;
Coherence_shuf.pval_coherence   = pval_coherence;

% z-score
info_obs        = info_perspike;
info_shuf       = SpInfo_shuf.info_shuf;
info_obs_persec = info_persec;
info_persec_shuf   = SpInfo_shuf.info_sec_shuf;

SpInfo_z = (info_obs - mean(info_shuf, 'omitmissing'))./std(info_shuf, 'omitmissing');
SpInfo_n = info_obs/mean(info_shuf, 'omitmissing');
SpInfo_sec_z = (info_obs_persec - mean(info_persec_shuf, 'omitmissing'))./std(info_persec_shuf, 'omitmissing');
SpInfo_sec_n = info_obs_persec/ mean(info_persec_shuf, 'omitmissing');

SpInfo_shuf.SpInfo_z        = SpInfo_z;
SpInfo_shuf.SpInfo_n        = SpInfo_n;
SpInfo_shuf.SpInfo_sec_z    = SpInfo_sec_z;
SpInfo_shuf.SpInfo_sec_n    = SpInfo_sec_n;


Cohe_Z = Coherence.Z;
% Coherence.SpCorr = Cohe_SpCorr;
Cohe_Z_shuf = Coherence_shuf.Z_shuf;
Cohe_Z_zscore = (Cohe_Z - mean(Cohe_Z_shuf, 'omitmissing'))./std(Cohe_Z_shuf, 'omitmissing');
Cohe_Z_n = Cohe_Z / mean(Cohe_Z_shuf, 'omitmissing');
Coherence_shuf.Z_zscore = Cohe_Z_zscore;
Coherence_shuf.Z_n      = Cohe_Z_n;

end

%% Spatial Coherence
function [Z, SpCorr, SmoothRateMap] = func_SpatialCoherence_251014(rateMap)
% https://www.jneurosci.org/content/9/12/4101
% https://pmc.ncbi.nlm.nih.gov/articles/PMC6282778/
%%
% Compute spatial coherence of a rate map (Muller 1989 method)
% INPUT:
%   rateMap : 2D matrix of firing rates [Hz]
%
% OUTPUT:
%   SpCorr         : Pearson correlation between each bin and its neighbors
%   Z              : Fisher's Z-transform of the correlation
%   SmoothRateMap  : 3×3 local mean map (center pixel excluded from mean)

if isempty(rateMap) || all(isnan(rateMap(:)))
    Z = NaN;
    SpCorr = NaN;
    SmoothRateMap = rateMap;
    return;
end

rateMap(isnan(rateMap)) = 0;

% Smooth with 3x3 boxcar
kernel = ones(5);
kernel(3,3) = 0;       % exclude center
kernel = kernel / sum(kernel(:));  % normalize (1/8)
SmoothRateMap = conv2(rateMap, kernel, 'same');

% Compute Pearson correlation
valid = ~isnan(rateMap) & ~isnan(SmoothRateMap);
if sum(valid(:)) < 5
    SpCorr = NaN;
    Z = NaN;
    return;
end
SpCorr = corr(rateMap(valid), SmoothRateMap(valid), 'Rows','complete');

% Fisher transform
Z = atanh(SpCorr);  % 0.5*log((1+R)/(1-R))

% figure
% imagesc(SmoothRateMap_forCoherence)
end

%% speed modulation
function speedModResult = func_analyzeSpeedModulation_251014(spikeFrames, smoothedSpeed, VideoSampling, doShuffle, nShuffle, doFitlm, doFitLogistic, minBinTime)
%%
% spikeFrames   : [Ns x 1] indeces of spike frames
% smoothedSpeed : [Nframes x 1] animal velocity at each frame [cm/s]
% VideoSampling : video sampling rate [Hz] (30Hz)
% Option:
%   'Shuffling' : 'on' / 'off' (default 'off')
%   'nShuffle'  : number of shuffling (default 1000)
%

%
frameRate = VideoSampling;
% ==== get spike frames and speed ====
spikeFrames = spikeFrames(~isnan(spikeFrames));
spikeSpeed = smoothedSpeed(spikeFrames);

% ==== firing rate vs. speed ====
edges = 0:20;
nBins = length(edges)-1;

SpeedCounts = histcounts(smoothedSpeed, edges);                % number of bin appearance
SpeedTimeHist   = SpeedCounts / frameRate;                           % convert from count to time
SpeedTimeHist(SpeedTimeHist==0) = NaN;

spikeCounts = histcounts(spikeSpeed, edges);
ratePerBin  = spikeCounts ./ SpeedTimeHist;                        % spike count hist to hist of firing rate [Hz]
speedCenters = edges(1:end-1) + diff(edges)/2;

% firing rate ==0 means that data is missing
ratePerBin(ratePerBin==0) = NaN;
ratePerBin_original = ratePerBin;
ratePerBin(SpeedTimeHist<minBinTime) = NaN;

% correlation
validCorr = sum(~isnan(ratePerBin)) >= 3;
if validCorr
    [corr_r, corr_p] = corr(speedCenters', ratePerBin', 'Rows','complete');
    [corr_SpearmanR, corr_SpearmanP] = corr(speedCenters', ratePerBin', 'Type','Spearman','Rows','complete');
else
    corr_r = NaN; corr_p = NaN;
    corr_SpearmanR = NaN; corr_SpearmanP = NaN;
end

%  linear regression
if doFitlm
    mdl = fitlm(speedCenters, ratePerBin);
    speedModResult.regression   = mdl;
    speedModResult.regression_basics.slope = mdl.Coefficients.Estimate(2);
    speedModResult.regression_basics.intercept = mdl.Coefficients.Estimate(1); % 切片
    speedModResult.regression_basics.p_slope = mdl.Coefficients.pValue(2); % whether slope > 0
    speedModResult.regression_basics.R2 = mdl.Rsquared.Ordinary;
    speedModResult.regression_basics.R2_adj = mdl.Rsquared.Adjusted;
    speedModResult.regression_basics.p_model = mdl.ModelFitVsNullModel.Pvalue;
    speedModResult.regression_basics.rmse = mdl.RMSE;
end

% sigmoid fit (Logistic regression)
if doFitLogistic
    fitLogistic = @(v,r) fit(v(:), r(:), ...
        fittype('L/(1+exp(-k*(x-x0)))','independent','x','coefficients',{'L','k','x0'}), ...
        'StartPoint',[max(r) 1 mean(v)]);
    v = speedCenters;
    r = ratePerBin;
    validIdx = ~isnan(v) & ~isnan(r);
    v_valid = v(validIdx);
    r_valid = r(validIdx);
    if numel(v_valid) >= 3  % calc only if there are at least three valid data points
        try
            f = fitLogistic(v_valid, r_valid);
            slopes = f.k;
        catch
            slopes = NaN;
        end
    else
        slopes = NaN;
    end
    speedModResult.LogisticSlopes = slopes;
end

%
speedModResult.speedInfo = func_calcSpeedInformation(ratePerBin, SpeedTimeHist);

% ==== shuffling (option) ====
shuffleR = [];
rateShuf = []; %#ok<NASGU>
p_shuffle = NaN;
zscoredRate = [];

if doShuffle && nShuffle>0
    nFrames = numel(smoothedSpeed);
    shuffleR = nan(nShuffle,1);
    shuffle_corr_SpearmanR = nan(nShuffle,1);
    rateShuf = nan(nBins, nShuffle);
    I_secShuf = nan(nShuffle,1);
    I_spikeShuf = nan(nShuffle,1);
    % tic
    % parfor i_shuffle = 1:nShuffle
    for i_shuffle = 1:nShuffle
        shift = randi(nFrames);
        spkShift = mod(spikeFrames + shift - 1, nFrames) + 1;
        local_smoothedSpeed = smoothedSpeed;
        spdShift = local_smoothedSpeed(spkShift);
        sc = histcounts(spdShift, edges);
        rpb = sc ./ SpeedTimeHist;
        shuffleR(i_shuffle) = corr(speedCenters', rpb', 'Rows','complete');
        shuffle_corr_SpearmanR(i_shuffle) = corr(speedCenters', rpb', 'Type', 'Spearman', 'Rows','complete');
        rateShuf(:,i_shuffle) = rpb;

        infoShuf = func_calcSpeedInformation(rpb, SpeedTimeHist);
        I_secShuf(i_shuffle) = infoShuf.I_sec;
        I_spikeShuf(i_shuffle) = infoShuf.I_spike;
    end
    % toc
    p_shuffle = (sum(abs(shuffleR) >= abs(corr_r)) + 1) / (nShuffle + 1);
    p_shuffle_Spearman = (sum(abs(shuffle_corr_SpearmanR) >= abs(corr_SpearmanR)) + 1) / (nShuffle + 1);
    shufMean = mean(rateShuf, 2, 'omitnan');
    shufStd  = std(rateShuf, 0, 2, 'omitnan');
    shufStd(shufStd==0) = NaN;
    zscoredRate    = (ratePerBin(:) - shufMean) ./ shufStd;

    speedModResult.p_shuffle_Spearman = p_shuffle_Spearman;
    speedModResult.shuffle_corr_SpearmanR = shuffle_corr_SpearmanR;

    I_secShuf_p = (sum(I_secShuf >= speedModResult.speedInfo.I_sec) + 1) / (nShuffle + 1);
    I_spikeShuf_p = (sum(I_spikeShuf >= speedModResult.speedInfo.I_spike) + 1) / (nShuffle + 1);


    % z-score
    info_obs        = speedModResult.speedInfo.I_spike;
    info_shuf       = I_spikeShuf;
    info_obs_persec = speedModResult.speedInfo.I_sec;
    info_persec_shuf   = I_secShuf;

    SpInfo_z = (info_obs - mean(info_shuf, 'omitmissing'))./std(info_shuf, 'omitmissing');
    SpInfo_n = info_obs/mean(info_shuf, 'omitmissing');
    SpInfo_sec_z = (info_obs_persec - mean(info_persec_shuf, 'omitmissing'))./std(info_persec_shuf, 'omitmissing');
    SpInfo_sec_n = info_obs_persec/ mean(info_persec_shuf, 'omitmissing');

    if std(info_shuf)==0, SpInfo_z = NaN; end
    if std(info_persec_shuf)==0, SpInfo_sec_z = NaN; end

    speedModResult.speedInfo.SpInfo_z        = SpInfo_z;
    speedModResult.speedInfo.SpInfo_n        = SpInfo_n;
    speedModResult.speedInfo.SpInfo_sec_z    = SpInfo_sec_z;
    speedModResult.speedInfo.SpInfo_sec_n    = SpInfo_sec_n;

    speedModResult.speedInfo.I_secShuf        = I_secShuf;
    speedModResult.speedInfo.I_secShuf_p       = I_secShuf_p;
    speedModResult.speedInfo.I_spikeShuf        = I_spikeShuf;
    speedModResult.speedInfo.I_spikeShuf_p        = I_spikeShuf_p;


end

% ==== summarize outputs ====
speedModResult.speedCenters = speedCenters;
speedModResult.ratePerBin   = ratePerBin;
speedModResult.ratePerBin_original = ratePerBin_original;
speedModResult.corr_r       = corr_r;
speedModResult.corr_p       = corr_p;
speedModResult.shuffle_r    = shuffleR;
speedModResult.shuffle_p    = p_shuffle;
speedModResult.zscoredRate  = zscoredRate;

speedModResult.corr_SpearmanR = corr_SpearmanR;
speedModResult.corr_SpearmanP = corr_SpearmanP;


end

% speed information
function speedInfo = func_calcSpeedInformation(ratePerBin, SpeedTimeHist)
% https://www.nature.com/articles/s41593-017-0055-3#Sec2
% calculate speed information (Skaggs-style)
% INPUT:
%   ratePerBin     : firing rate per speed bin [Hz]
%   SpeedTimeHist  : time spent in each speed bin [s]
% OUTPUT:
%   speedInfo.I_sec    : bits/sec
%   speedInfo.I_spike  : bits/spike
%   speedInfo.meanRate : mean firing rate [Hz]

validBins = SpeedTimeHist > 0;
p_v = SpeedTimeHist(validBins) / sum(SpeedTimeHist(validBins)); % occupancy prob
lambda_v = ratePerBin(validBins);

meanRate = sum(p_v .* lambda_v, 'omitmissing');

% avoid log2(0)
nonzero = lambda_v > 0 & meanRate > 0;
term = zeros(size(lambda_v));
term(nonzero) = p_v(nonzero) .* lambda_v(nonzero) .* log2(lambda_v(nonzero)/meanRate);

I_sec = sum(term, 'omitmissing');
I_spike = I_sec / meanRate;

speedInfo.I_sec    = I_sec;
speedInfo.I_spike  = I_spike;
speedInfo.meanRate = meanRate;
end


%% stability and reliability
function [fig2, HalveCompareField, HalveCompareSpeed] = ...
    func_HalveCompare(TraceData, allFrame_vid, Time_aligned, mousePos, spikeFrames, spikePos, SpeedResults,...
    binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, nThresholds_speed, CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, speedInd_list,...
    doImShow, doCompare)
% Within-trial stability was the correlation between the spatially corresponding bins of the rate maps from the temporal first and last halves of the trial.
% https://www.cell.com/neuron/fulltext/S0896-6273(15)00418-3
% https://pmc.ncbi.nlm.nih.gov/articles/PMC3547638/#sec2
%%
% doImShow = 1;
% doCompare = 'Stability';
% doCompare = 'Reliability';
% speedInd_list = [1 2 3 4 5];

Fs = size(TraceData,1)/TraceData(end,1);

thrVals = SpeedResults.Thresholds(speedInd_list);
n_speedInd = length(speedInd_list);

occThresh_halve = 1;

rateMapUniform_sp_halve = cell(2,1);
rateMapUniform_Gaussian_sp_halve = cell(2,1);
spikeMap_Rawcount_sp_halve = cell(2,1);
numEvent_sp_halve = cell(2,1);
runSpkIdx_sp_halve = cell(2,1);
neuronPos_sp_halve = cell(2,1);
occMapUniform_Gaussian_sp_halve = cell(2,1);
occMapUniform_sp_halve = cell(2,1);
for i = 1:2
    rateMapUniform_sp_halve{i}          = cell(nThresholds_speed, 1);
    rateMapUniform_Gaussian_sp_halve{i}        = cell(nThresholds_speed, 1);
    spikeMap_Rawcount_sp_halve{i} = cell(nThresholds_speed, 1);
    numEvent_sp_halve{i}         = nan(nThresholds_speed, 1);
    runSpkIdx_sp_halve{i}        = cell(nThresholds_speed, 1);
    neuronPos_sp_halve{i}        = cell(nThresholds_speed, 1);
    occMapUniform_Gaussian_sp_halve{i}          = cell(nThresholds_speed, 1);
    occMapUniform_sp_halve{i}         = cell(nThresholds_speed, 1);
end

runFrames_halve = cell(2,1);
Cor_old = cell(2,n_speedInd);
Cor_new = cell(2,n_speedInd);
RunningTime = nan(n_speedInd, 2);
for thr_sp_ind = 1:n_speedInd
    thr_sp = speedInd_list(thr_sp_ind);

    % spike position during running
    runFrames = SpeedResults.idx{thr_sp}; %Index of running frame, extracted from allFrame_vid

    if strcmp(doCompare, 'Stability')
        % divide running frames by halves of all recording
        midFrame = round(length(allFrame_vid)/2);
        midTime = Time_aligned(midFrame);
        runFrames_halve{1}  = runFrames(runFrames <= midFrame);
        runFrames_halve{2} = runFrames(runFrames > midFrame);
    elseif strcmp(doCompare, 'Reliability')
        runFrames_halve{1}  = runFrames(mod(runFrames,2) == 0); %even frames
        runFrames_halve{2} = runFrames(mod(runFrames,2) == 1);  %odd frames
    end

    for i = 1:2
        runSpkIdx_halve = ismember(spikeFrames, runFrames_halve{i});
        neuronPos_halve  = spikePos(runSpkIdx_halve, :);
        mousePos_run_halve  = mousePos(runFrames_halve{i}, :);

        mousePos_run = mousePos_run_halve;
        neuronPos = neuronPos_halve;

        % calculate rate map using adaptive binning
        [occMapMerged, spkMapMerged, rateMapMerged, xEdges, yEdges, nBin, Dur_Running, areaMap, rateMapMerged_ori, occMapMerged_ori] = ...
            func_adaptiveBinRateMap(mousePos_run, neuronPos, AnimalTrack, binSize_cm, VideoSampling, occThresh_halve, sigma, kernelSize); %#ok<ASGLU>

        Edges = {xEdges, yEdges};
        [occMapNormalized, occMapGaussian, occMapUniform, occMapUniform_Gaussian, X_old, Y_old, X_new, Y_new, xPlot, yPlot, xDim, yDim] = ...
            computeGaussianOccMap(AnimalTrack, occMapMerged, Edges, areaMap, binSize_cm, sigma, kernelSize); %#ok<ASGLU>
        [rateMapUniform, rateMapUniform_Gaussian] = computeGaussianRateMap(rateMapMerged, X_old, Y_old, X_new, Y_new, sigma, kernelSize);

        Cor_old{i, thr_sp_ind} = {X_old, Y_old};
        Cor_new{i, thr_sp_ind} = {X_new, Y_new};

        % --- save variables ---
        rateMapUniform_sp_halve{i}{thr_sp_ind, 1}          = rateMapUniform;
        rateMapUniform_Gaussian_sp_halve{i}{thr_sp_ind, 1}         = rateMapUniform_Gaussian;
        spikeMap_Rawcount_sp_halve{i}{thr_sp_ind, 1} = spkMapMerged;
        numEvent_sp_halve{i}(thr_sp_ind, 1)         = sum(spkMapMerged(:),'omitmissing');
        runSpkIdx_sp_halve{i}{thr_sp_ind, 1}        = runSpkIdx_halve;
        neuronPos_sp_halve{i}{thr_sp_ind, 1}        = neuronPos_halve;
        occMapUniform_Gaussian_sp_halve{i}{thr_sp_ind, 1}          = occMapUniform_Gaussian;
        occMapUniform_sp_halve{i}{thr_sp_ind, 1}         = occMapUniform;
        RunningTime(thr_sp_ind, i) = length(runFrames_halve{i})/Fs;

        % figure
        % pcolor(X_old, Y_old, spkMapMerged)
    end
end
% toc %0.086sec


fig2 = [];
if doImShow
    % fig2 = figure('Position',[-1897          72        1883         876]);
    fig2 = figure('Position',[-1897          72        1883         876],'Visible','off');
    % get(gcf,'Position')
end

% arena size
cropRect = AnimalTrack.TrackData.cropRect;
psh = AnimalTrack.TrackData.psh; psw = AnimalTrack.TrackData.psw;
xDim = cropRect(3)*psw; yDim = cropRect(4)*psh;
xmin = 0; xmax = xDim; ymin = 0; ymax = yDim;

% Firing rates of all animal speeds
TraceData_halve = cell(2,1);
RecordingDur = nan(2,1);
FiringRate_halve = nan(2,1);
if strcmp(doCompare, 'Stability')
    [~, idxClosest] = min(abs(TraceData(:,1) - midTime));
    Frames_CaTrace_halve{1} = 1:idxClosest;
    Frames_CaTrace_halve{2} = idxClosest:size(TraceData,1);
    for i=1:2
        TraceData_halve{i} = TraceData(Frames_CaTrace_halve{i},:);
        [NonRecordingDur, ~] = func_detectNonRecording(TraceData_halve{i}(:,2), Fs);
        RecordingDur(i) = (TraceData_halve{i}(end,1) - TraceData_halve{i}(1,1) - NonRecordingDur);
        FiringRate_halve(i) = nnz(runSpkIdx_sp_halve{i}{1, 1})/RecordingDur(i);
    end
elseif strcmp(doCompare, 'Reliability')
    for i=1:2
        FiringRate_halve(i) = nnz(runSpkIdx_sp_halve{i}{1, 1})/(TraceData(end,1)/2); % it does not count for non-recording period
    end
end

% --- Show Calcium trace and Occupancy maps ---
if doImShow
    if strcmp(doCompare, 'Reliability')
        midTime = [];
    end
    func_showCatraceAndOccMap_forStability(TraceData, midTime, CaProcessingParm, occMapUniform_Gaussian_sp_halve, runSpkIdx_sp_halve, i_Cell, thresholds_speed, FiringRate_halve, doCompare)
end

%preallocation
speedModResult_halve = cell(2,1);
meanRate_info = cell(2,n_speedInd);
info_spike = cell(2,n_speedInd);
info_sec = cell(2,n_speedInd);
info_spatialSparsity = cell(2,n_speedInd);
Coherence_Z = cell(2,n_speedInd);
Coherence_SpCorr = cell(2,n_speedInd);
for i=1:2
    if strcmp(doCompare, 'Stability')
        if i == 1
            Frames_Video_halve = 1:midFrame;
        elseif i == 2
            Frames_Video_halve = midFrame:length(allFrame_vid);
        end
    elseif strcmp(doCompare, 'Reliability')
        if i == 1
            Frames_Video_halve  = allFrame_vid(mod(allFrame_vid,2) == 0); %even frames
        elseif i == 2
            Frames_Video_halve  = allFrame_vid(mod(allFrame_vid,2) == 1); %odd frames
        end
    end

    x_animalall_halve = mousePos(Frames_Video_halve,1);
    y_animakall_halve = mousePos(Frames_Video_halve,2);

    % --- Show AnimalTrace and SpikePos ---
    if doImShow
        func_showAnimalTraceAndSpikePos_stability_251015(speedInd_list, runSpkIdx_sp_halve, i, neuronPos_sp_halve, x_animalall_halve, y_animakall_halve, xmin,xmax,ymin,ymax,thrVals, RunningTime)
    end
end


%% calc spatial info
rateMapUniform_Gaussian_sp_halve_all = cell(2,n_speedInd);
rateMapUniform_sp_halve_all = cell(2,n_speedInd);
occMapUniform_Gaussian_sp_halve_all = cell(2,n_speedInd);
for i=1:2
    PeakFiringRate_sp = cell(n_speedInd,1);
    for k = 1:n_speedInd

        rateMapUniform_Gaussian = rateMapUniform_Gaussian_sp_halve{i}{k, 1};
        rateMapUniform = rateMapUniform_sp_halve{i}{k, 1};

        if doImShow
            subplot(5,10,5+k + 10*(i - 1)),
            pcolor(X_new, Y_new, rateMapUniform); shading flat;
            axis tight equal off; set(gca, 'YDir','reverse')
            PeakFR = max(rateMapUniform,[],'all','omitnan');
            if i == 1
                title({'Early raw rate map', sprintf('Peak firing rate = %.3f Hz', PeakFR)})
            else
                title({'Later raw rate map', sprintf('Peak firing rate = %.3f Hz', PeakFR)})
            end

            subplot(5,10,25+k + 10*(i - 1))
            pcolor(X_new, Y_new, rateMapUniform_Gaussian); shading flat;
            axis tight equal off; set(gca, 'YDir','reverse')
            PeakFiringRate_sp{k} = max(rateMapUniform_Gaussian,[],'all','omitnan');
            if i == 1
                title({'Early gaussian rate map', sprintf('Peak firing rate = %.3f Hz', PeakFiringRate_sp{k})})
            else
                title({'Later gaussian rate map', sprintf('Peak firing rate = %.3f Hz', PeakFiringRate_sp{k})})
            end

        end

        occ_temp = occMapUniform_Gaussian_sp_halve{i}{k, 1};
        [meanRate_info{i,k}, info_spike{i,k}, info_sec{i,k}, info_spatialSparsity{i,k}] = func_SkaggsSpatialInfo_251014(rateMapUniform_Gaussian, occ_temp);
        [Coherence_Z{i,k}, Coherence_SpCorr{i,k}, ~] = func_SpatialCoherence_251014(rateMapUniform);

        rateMapUniform_Gaussian_sp_halve_all{i,k}   = rateMapUniform_Gaussian;
        rateMapUniform_sp_halve_all{i,k}            = rateMapUniform;
        occMapUniform_Gaussian_sp_halve_all{i,k}    = occ_temp;
    end

end

%
HalveCompareField.SpInfo.meanRate = meanRate_info;
HalveCompareField.SpInfo.info_spike = info_spike;
HalveCompareField.SpInfo.info_sec = info_sec;
HalveCompareField.SpInfo.info_spatialSparsity = info_spatialSparsity;
HalveCompareField.SpInfo.Coherence_Z = Coherence_Z;
HalveCompareField.SpInfo.Coherence_SpCorr = Coherence_SpCorr;


%% calc stability (Pearson R)

% calc. stability parameters
HalvesCorr_Pearson_sp = nan(n_speedInd,1);
HalvesCorr_Pearson_p_sp = nan(n_speedInd,1);
HalvesCorr_Pearson_z_sp = nan(n_speedInd,1);
HalvesCorr_Spearman_sp = nan(n_speedInd,1);
HalvesCorr_Spearman_p_sp = nan(n_speedInd,1);
XCor_sp = cell(n_speedInd,1);
XCor_max_sp = nan(n_speedInd,1);
shift_x_cm_sp = nan(n_speedInd,1);
shift_y_cm_sp = nan(n_speedInd,1);
shiftDist_sp = nan(n_speedInd,1);
for k = 1:n_speedInd
    A = rateMapUniform_Gaussian_sp_halve{1}{k, 1};
    B = rateMapUniform_Gaussian_sp_halve{2}{k, 1};
    A(isnan(A))=0; B(isnan(B))=0;

    % stability/reliability based on correlation
    [HalvesCorr_Pearson, p_Pearson] = corr(A(:), B(:));
    r = HalvesCorr_Pearson;
    z_Pearson = 0.5 * log((1 + r) / (1 - r)); %Fisher z ( = atanh(r))
    [HalvesCorr_Spearman, p_Spearman] = corr(A(:), B(:), 'Type', 'Spearman');

    % stability based on cross-correlation shift
    if std(A(:)) == 0 || std(B(:)) == 0
        % flat template or image → set outputs to zero
        shift_x_cm = NaN;
        shift_y_cm = NaN;
        shiftDist = NaN;
        XCor = zeros(size(A)+size(B)-1); % normxcorr2 would return this size

        dx = 1; dy = 1;
        XCorPeak = NaN;
    else
        % normal cross-correlation
        XCor = normxcorr2(A, B);
        [XCorPeak, imax] = max(XCor(:));
        [ypeak, xpeak] = ind2sub(size(XCor), imax);
        off_y = ypeak - size(A,1);
        off_x = xpeak - size(A,2);
        dx =  binSize_cm/factor;
        dy =  binSize_cm/factor;
        shift_x_cm = off_x * dx;
        shift_y_cm = off_y * dy;
        shiftDist = hypot(shift_x_cm, shift_y_cm);
    end

    if doImShow
        subplot(5,10,45+k)
        ny = size(XCor,1);  nx = size(XCor,2);
        yAxis = (1:ny) - size(A,1);  % shift toward row direction
        xAxis = (1:nx) - size(A,2);  % shift toward colomn direction
        xAxis_cm = xAxis * dx; yAxis_cm = yAxis * dy;

        imagesc(xAxis_cm, yAxis_cm, XCor);
        axis equal off, colormap jet
        axis off equal tight
        axis xy; % zero is bottom in the y axis
        s_corr = sprintf("X-Cor map; Cor = %.2f", HalvesCorr_Pearson);
        s_xcorr = sprintf("XCor = %.2f, Shift = %.2fcm", XCorPeak, shiftDist);
        title({s_corr, s_xcorr});
        clim([0 0.8])
    end

    HalvesCorr_Pearson_sp(k) = HalvesCorr_Pearson;
    HalvesCorr_Pearson_p_sp(k) = p_Pearson;
    HalvesCorr_Pearson_z_sp(k) = z_Pearson;

    HalvesCorr_Spearman_sp(k) = HalvesCorr_Spearman;
    HalvesCorr_Spearman_p_sp(k) = p_Spearman;

    XCor_sp{k} = XCor;
    XCor_max_sp(k) = XCorPeak;
    shift_x_cm_sp(k) = shift_x_cm;
    shift_y_cm_sp(k) = shift_y_cm;
    shiftDist_sp(k) = shiftDist;
end

% summarize data
HalveCompareField.speedInd_list = speedInd_list;

HalveCompareField.RunningTime = RunningTime;
HalveCompareField.numEvent_sp_halve = numEvent_sp_halve;
HalveCompareField.FiringRate_halve = FiringRate_halve;

HalveCompareField.rateMapUniform_Gaussian_sp_halve_all  = rateMapUniform_Gaussian_sp_halve_all;
HalveCompareField.rateMapUniform_sp_halve_all           = rateMapUniform_sp_halve_all;
HalveCompareField.occMapUniform_Gaussian_sp_halve_all   = occMapUniform_Gaussian_sp_halve_all;

HalveCompareField.HalvesCorr_Pearson_sp = HalvesCorr_Pearson_sp;
HalveCompareField.HalvesCorr_Pearson_p_sp = HalvesCorr_Pearson_p_sp;
HalveCompareField.HalvesCorr_Pearson_z_sp = HalvesCorr_Pearson_z_sp;
HalveCompareField.HalvesCorr_Spearman_sp = HalvesCorr_Spearman_sp;
HalveCompareField.HalvesCorr_Spearman_p_sp = HalvesCorr_Spearman_p_sp;

HalveCompareField.XCor_sp = XCor_sp;
HalveCompareField.XCor_max_sp = XCor_max_sp;
HalveCompareField.XCor_shift_x_cm_sp = shift_x_cm_sp;
HalveCompareField.XCor_shift_y_cm_sp = shift_y_cm_sp;
HalveCompareField.XCor_shiftDist_sp = shiftDist_sp;


%% calc speed modulation stability/reliablity
ratePerBin_halves = nan(2,20);

% %sigmoid fit
% fitLogistic = @(v,r) fit(v(:), r(:), ...
%     fittype('L/(1+exp(-k*(x-x0)))','independent','x','coefficients',{'L','k','x0'}), ...
%     'StartPoint',[max(r) 1 mean(v)]);

slopes = nan(1,2);
for i = 1:2
    if strcmp(doCompare, 'Stability')
        if i == 1
            spikeFrames_halve  = spikeFrames(spikeFrames <= midFrame);
            smoothedSpeed_halve = smoothedSpeed;
            smoothedSpeed_halve(midFrame:end) = NaN;
        elseif i == 2
            spikeFrames_halve = spikeFrames(spikeFrames > midFrame);
            smoothedSpeed_halve = smoothedSpeed;
            smoothedSpeed_halve(1:midFrame) = NaN;
        end
    elseif strcmp(doCompare, 'Reliability')
        if i == 1
            spikeFrames_halve  = spikeFrames(mod(spikeFrames,2) == 0); %even frames
            smoothedSpeed_halve = smoothedSpeed;
            % smoothedSpeed_halve(mod(spikeFrames,2) == 1) = NaN;
            smoothedSpeed_halve(mod(allFrame_vid,2) == 1) = NaN;
        elseif i == 2
            spikeFrames_halve  = spikeFrames(mod(spikeFrames,2) == 1); %odd frames
            smoothedSpeed_halve = smoothedSpeed;
            % smoothedSpeed_halve(mod(spikeFrames,2) == 0) = NaN;
            smoothedSpeed_halve(mod(allFrame_vid,2) == 0) = NaN;
        end

    end
    doShuffle = 0; nShuffle = 0; doFitlm = 0; doFitLogistic = 1;
    minBinTime = 30;
    speedModResult_halve{i} = func_analyzeSpeedModulation_251014(spikeFrames_halve, smoothedSpeed_halve, VideoSampling, doShuffle, nShuffle, doFitlm, doFitLogistic, minBinTime);
    ratePerBin = speedModResult_halve{i}.ratePerBin;
    ratePerBin_halves(i,1:20) = ratePerBin;
    speedCenters = speedModResult_halve{i}.speedCenters;

    slopes(i) = speedModResult_halve{i}.LogisticSlopes;
end
[r_Pearson, p_Pearson] = corr(ratePerBin_halves(1,:)', ratePerBin_halves(2,:)', 'Rows','complete');
[r_Spearman, p_Spearman] = corr(ratePerBin_halves(1,:)', ratePerBin_halves(2,:)', 'Rows','complete', 'Type', 'Spearman');
stability_logistic = 1 - abs(diff(slopes))/max(abs(slopes));  % 0: unstable, 1:stable

if doImShow
    subplot(5,10,5)
    b = bar(speedCenters, ratePerBin_halves);
    for i = 1:2
        b(i).EdgeColor = 'none';
        b(i).BarWidth = 1.5;
        b(i).FaceAlpha = 0.7;
    end
    xlabel('Speed (cm/s)');
    ylabel('Firing rate (Hz)');
    title({'Speed modulation stability', ...
        sprintf( 'Pearson r = %.2f, p = %.3f', r_Pearson, p_Pearson), sprintf( 'Spearman r = %.2f, p = %.3f', r_Spearman, p_Spearman), ...
        sprintf("slope diff. of logistic = %.2f", stability_logistic)});
    xlim([0 15])
    box off
    hold on
    if strcmp(doCompare, 'Stability')
        sgtitle(strcat(s_figtitle, " stability"));
    elseif strcmp(doCompare, 'Reliability')
        sgtitle(strcat(s_figtitle, " reliability"));
    end
end

% summarize data
HalveCompareSpeed.ratePerBin_halves = ratePerBin_halves;
HalveCompareSpeed.speedCenters = speedCenters;
HalveCompareSpeed.r_Pearson = r_Pearson;
HalveCompareSpeed.p_Pearson = p_Pearson;
HalveCompareSpeed.r_Spearman = r_Spearman;
HalveCompareSpeed.p_Spearman = p_Spearman;
HalveCompareSpeed.stability_logistic = stability_logistic;

end

% ---- Helper functions ----
function func_showCatraceAndOccMap_forStability(TraceData, midTime, CaProcessingParm, occMapSmooth_sp_halve, runSpkIdx_sp_halve, i_Cell, thresholds_speed, FiringRate_halve, doCompare)

TraceData_halve = cell(2,1);
for i=1:2

    if strcmp(doCompare, 'Stability')
        [~, idxClosest] = min(abs(TraceData(:,1) - midTime));
        Frames_CaTrace_halve{1} = 1:idxClosest;
        Frames_CaTrace_halve{2} = idxClosest:size(TraceData,1);

        TraceData_halve{i} = TraceData(Frames_CaTrace_halve{i},:);

        subplot(5,10,[1 2] + 2*(i-1))
        plot(TraceData_halve{i}(:,1), TraceData_halve{i}(:,i_Cell+1),'k'), axis tight, box off
        hold on, yline(CaProcessingParm.threshold(i_Cell),'--','Color',[1 0 0 0.6])
        xlabel('time [s]')
        if i == 1
            title(sprintf('Earlier Firing rate = %.3f Hz; %d events', FiringRate_halve(i), nnz(runSpkIdx_sp_halve{i}{1, 1})))
        elseif i == 2
            title(sprintf('Later Firing rate = %.3f Hz; %d events', FiringRate_halve(i), nnz(runSpkIdx_sp_halve{i}{1, 1})))
        end
    end

    % Show occupancy map
    for k = 1:5
        subplot(5,10,30 + k + 10*(i - 1))
        occMapSmooth = occMapSmooth_sp_halve{i}{k, 1};
        imagesc(occMapSmooth)
        axis tight equal off
        if i == 1
            title(sprintf('Early occ map at %dcm/s', thresholds_speed(k)))
        else
            title(sprintf('Later occ map at %dcm/s', thresholds_speed(k)))
        end
        colormap jet
        if max(occMapSmooth(:))>0
            clim([0 max(occMapSmooth(:))])
        end
    end

end

end

function func_showAnimalTraceAndSpikePos_stability_251015(speedInd_list,runSpkIdx_sp_halve, i, neuronPos_sp_halve, x_animalall_halve, y_animakall_halve, xmin,xmax,ymin,ymax,thrVals, RunningTime)
for k = 1:length(speedInd_list) %1:5
    spId   = speedInd_list(k);
    runIdx = runSpkIdx_sp_halve{i}{spId, 1};
    % spike scatter
    subplot(5,10, 10 + k + 10*(i - 1))
    x_spk = neuronPos_sp_halve{i}{k, 1}(:,1);
    y_spk = neuronPos_sp_halve{i}{k, 1}(:,2);
    plotScatter(x_animalall_halve, y_animakall_halve, x_spk, y_spk, xmin,xmax,ymin,ymax,thrVals(k),nnz(runIdx))
    if i == 1
        title({sprintf('Early >%.1fcm/s, %d events', thrVals(k),nnz(runIdx)), sprintf('Run time: %.1f min', RunningTime(k,i)/60)})
    else
        title({sprintf('Later >%.1fcm/s, %d events', thrVals(k),nnz(runIdx)), sprintf('Run time: %.1f min', RunningTime(k,i)/60)})
    end
end
end



%% Shuffling stability and reliability

function Halves_shuf = func_StabilityShuffle(TraceData, allFrame_vid, Time_aligned, mousePos, ...
    spikeFrames, spikePos, SpeedResults, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, ...
    nThresholds_speed, CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, doCompare, shuffIdx)
%% initialize parameters
speedInd_list = 3;
doImShow = 0; % 0: DO NOT SHOW IMAGE
% doCompare = 'Stability';
nShuffle = 1000; % change here

occThresh_halve = 1;


%% first, compute real stability or reliability
[~, Field_real, Speed_real] = func_HalveCompare( ...
    TraceData, allFrame_vid, Time_aligned, mousePos, spikeFrames, spikePos, SpeedResults, ...
    binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, nThresholds_speed, ...
    CaProcessingParm, i_Cell, thresholds_speed, factor, smoothedSpeed, s_figtitle, ...
    speedInd_list, doImShow, doCompare);

Halves_shuf.Field.real = Field_real;
Halves_shuf.Speed.real = Speed_real;


%% start shuffling
% close all

thr_sp = speedInd_list; % use this for temporarily, this might be changed in future
runFrames = SpeedResults.idx{thr_sp};

mousePosTime = [(1:length(mousePos(:,1)))', Time_aligned', mousePos];
mousePosTime_run = mousePosTime(runFrames,:);

midFrame = round(length(mousePosTime(:,1))/2);
midTime = Time_aligned(midFrame);

%Gaussian Occupancy Map
if strcmp(doCompare, 'Stability')
    runFrames_halve{1}  = runFrames(runFrames <= midFrame);
    runFrames_halve{2} = runFrames(runFrames > midFrame);
elseif strcmp(doCompare, 'Reliability')
    runFrames_halve{1}  = runFrames(mod(runFrames,2) == 0); %even frames
    runFrames_halve{2} = runFrames(mod(runFrames,2) == 1);  %odd frames
end

occMapUniform_Gaussian_halve_ForShuf = cell(2,1);
mousePos_run_halve_ForShuf = cell(2,1);
Edges_ForShuf  = cell(2,1);
XY_old  = cell(2,1);
XY_new  = cell(2,1);
for i = 1:2
    mousePos_run_halve  = mousePos(runFrames_halve{i}, :);

    neuronPos_dummy = [1,1];
    [occMapMerged, ~, ~, xEdges, yEdges, ~, ~, areaMap, ~, ~] = ...
        func_adaptiveBinRateMap(mousePos_run_halve, neuronPos_dummy, AnimalTrack, binSize_cm, VideoSampling, occThresh_halve, sigma, kernelSize);
    Edges_ForShuf{i}  = {xEdges, yEdges};
    [~, ~, ~, occMapUniform_Gaussian, X_old, Y_old, X_new, Y_new, ~, ~, ~, ~] = ...
        computeGaussianOccMap(AnimalTrack, occMapMerged, Edges_ForShuf{i}, areaMap, binSize_cm, sigma, kernelSize);

    XY_old{i} = {X_old, Y_old};
    XY_new{i} = {X_new, Y_new};


    mousePos_run_halve_ForShuf{i} = mousePos_run_halve;
    occMapUniform_Gaussian_halve_ForShuf{i} = occMapUniform_Gaussian;

end


%% shuffling from here
% preallocation
HalvesCorr_Pearson_shuf = nan(nShuffle, 1);
HalvesCorr_Spearman_shuf = nan(nShuffle, 1);
HalvesXCorMax_shuf = nan(nShuffle, 1);
HalvesXCorshiftDist_shuf = nan(nShuffle, 1);

meanRate_shuf = nan(nShuffle, 2);
info_shuf = nan(nShuffle, 2);
info_sec_shuf = nan(nShuffle, 2);
spatialSparsity_shuf = nan(nShuffle, 2);
Z_shuf = nan(nShuffle, 2);
Coherence_SpCorr = nan(nShuffle, 2);

% tic
% figure
% for i_shuffle = 1:nShuffle
parfor i_shuffle = 1:nShuffle

    ShufIdx = shuffIdx{i_shuffle};
    mousePosTime_run_local = mousePosTime_run;
    neuronPosTime_shuff_temp = mousePosTime_run_local(ShufIdx,:);
    neuronPosTime_shuff = cell(2,1);
    if strcmp(doCompare, 'Stability')
        [~, Ind] = min(abs(neuronPosTime_shuff_temp(:,2) - midTime));
        neuronPosTime_shuff{1} = neuronPosTime_shuff_temp(1:Ind-1, :);
        neuronPosTime_shuff{2} = neuronPosTime_shuff_temp(Ind:end, :);
    elseif strcmp(doCompare, 'Reliability')
        Ind = mod(neuronPosTime_shuff_temp(:,1),2) == 0;
        neuronPosTime_shuff{1} = neuronPosTime_shuff_temp(Ind, :); %even
        Ind = mod(neuronPosTime_shuff_temp(:,1),2) == 1;
        neuronPosTime_shuff{2} = neuronPosTime_shuff_temp(Ind, :); %odd
    end

    rateMapUniform_Gaussian_shuf_halve = cell(2,1);
    for i = 1:2
        % compute spike and rate map in shuffling loop
        neuronPos_shuff_temp = neuronPosTime_shuff{i}(:,3:4);

        local_mousePos_run_halve_ForShuf = mousePos_run_halve_ForShuf;
        mousePos_run = local_mousePos_run_halve_ForShuf{i};

        [~, ~, rMap_shuf, ~, ~, ~, ~, ~, ~, ~] = func_adaptiveBinRateMap(mousePos_run, neuronPos_shuff_temp, AnimalTrack, binSize_cm, VideoSampling, occThresh_halve, sigma, kernelSize);

        local_XY_old = XY_old;
        local_XY_new = XY_new;
        X_old =  local_XY_old{i}{1}; Y_old =  local_XY_old{i}{2};
        X_new =  local_XY_new{i}{1}; Y_new =  local_XY_new{i}{2};
        [rateMapUniform_shuf, rateMapUniform_Gaussian_shuf] = computeGaussianRateMap(rMap_shuf, X_old, Y_old, X_new, Y_new, sigma, kernelSize);

        local_occMapUniform_Gaussian_halve_ForShuf = occMapUniform_Gaussian_halve_ForShuf;
        occMapUniform_Gaussian = local_occMapUniform_Gaussian_halve_ForShuf{i};
        [meanRate_shuf(i_shuffle, i), info_shuf(i_shuffle, i), info_sec_shuf(i_shuffle, i), spatialSparsity_shuf(i_shuffle, i)] = func_SkaggsSpatialInfo_251014(rateMapUniform_Gaussian_shuf, occMapUniform_Gaussian);
        [Z_shuf(i_shuffle, i), Coherence_SpCorr(i_shuffle, i), ~] = func_SpatialCoherence_251014(rateMapUniform_shuf);

        rateMapUniform_Gaussian_shuf_halve{i, 1} = rateMapUniform_Gaussian_shuf;

        % subplot(1,2,i)
        % imagesc(rateMapUniform_Gaussian_shuf)
    end

    A = rateMapUniform_Gaussian_shuf_halve{1, 1};     B = rateMapUniform_Gaussian_shuf_halve{2, 1};
    A(isnan(A))=0; B(isnan(B))=0;

    % stability based on correlation
    [HalvesCorr_Pearson, ~] = corr(A(:), B(:));
    [HalvesCorr_Spearman, ~] = corr(A(:), B(:), 'Type', 'Spearman');

    % stability based on cross-correlation shift
    if std(A(:)) == 0 || std(B(:)) == 0
        % flat template or image → set outputs to zero
        shiftDist = 0;
        XCor = zeros(size(A)+size(B)-1); % normxcorr2 would return this size
    else
        % normal cross-correlation
        XCor = normxcorr2(A, B);
        [~, imax] = max(XCor(:));
        [ypeak, xpeak] = ind2sub(size(XCor), imax);
        off_y = ypeak - size(A,1);
        off_x = xpeak - size(A,2);
        dx =  binSize_cm/factor;
        dy =  binSize_cm/factor;
        shift_x_cm = off_x * dx;
        shift_y_cm = off_y * dy;
        shiftDist = hypot(shift_x_cm, shift_y_cm);
    end

    HalvesCorr_Pearson_shuf(i_shuffle) = HalvesCorr_Pearson;
    HalvesCorr_Spearman_shuf(i_shuffle) = HalvesCorr_Spearman;
    HalvesXCorMax_shuf(i_shuffle) = max(XCor(:));
    HalvesXCorshiftDist_shuf(i_shuffle) = shiftDist;

end
% toc

%%
% p test (one-side)
obs = cell2mat(Halves_shuf.Field.real.SpInfo.info_spike)';
shuf = info_shuf;
pval_spinfo     = (sum(shuf >= obs) + 1) / (nShuffle + 1);

obs = cell2mat(Halves_shuf.Field.real.SpInfo.info_sec)';
shuf = info_sec_shuf;
pval_spsecinfo   = (sum(shuf >= obs) + 1) / (nShuffle + 1);

obs = cell2mat(Halves_shuf.Field.real.SpInfo.info_spatialSparsity)';
shuf = spatialSparsity_shuf;
pval_sparsity   = (sum(shuf <= obs) + 1) / (nShuffle + 1);

obs = cell2mat(Halves_shuf.Field.real.SpInfo.Coherence_Z)';
shuf = Z_shuf;
pval_coherence   = (sum(shuf >= obs) + 1) / (nShuffle + 1);


%%
Halves_shuf.Field.shuf.meanRate_shuf = meanRate_shuf;
Halves_shuf.Field.shuf.info_shuf = info_shuf;
Halves_shuf.Field.shuf.info_sec_shuf = info_sec_shuf;
Halves_shuf.Field.shuf.spatialSparsity_shuf = spatialSparsity_shuf;
Halves_shuf.Field.shuf.Z_shuf = Z_shuf;
Halves_shuf.Field.shuf.Coherence_SpCorr = Coherence_SpCorr;

Halves_shuf.Field.shuf.pval_spinfo = pval_spinfo;
Halves_shuf.Field.shuf.pval_spsecinfo = pval_spsecinfo;
Halves_shuf.Field.shuf.pval_sparsity = pval_sparsity;
Halves_shuf.Field.shuf.pval_coherence = pval_coherence;

Halves_shuf.Field.shuf.HalvesCorr_Pearson_shuf      = HalvesCorr_Pearson_shuf;
Halves_shuf.Field.shuf.HalvesCorr_Spearman_shuf     = HalvesCorr_Spearman_shuf;
Halves_shuf.Field.shuf.HalvesXCorMax_shuf           = HalvesXCorMax_shuf;
Halves_shuf.Field.shuf.HalvesXCorshiftDist_shuf     = HalvesXCorshiftDist_shuf;

%%
% test (one-side)
obs     = Halves_shuf.Field.real.HalvesCorr_Pearson_sp;
shuf    = Halves_shuf.Field.shuf.HalvesCorr_Pearson_shuf;
Halves_shuf.Field.shuf.HalvesCorr_Pearson_shuf_pval    = (sum(shuf >= obs) + 1) / (nShuffle + 1);

obs     = Halves_shuf.Field.real.HalvesCorr_Spearman_sp;
shuf    = Halves_shuf.Field.shuf.HalvesCorr_Spearman_shuf ;
Halves_shuf.Field.shuf.HalvesCorr_Spearman_shuf_pval    = (sum(shuf >= obs) + 1) / (nShuffle + 1);

obs     = Halves_shuf.Field.real.XCor_max_sp;
shuf    = Halves_shuf.Field.shuf.HalvesXCorMax_shuf;
Halves_shuf.Field.shuf.HalvesXCorMax_shuf_pval    = (sum(shuf >= obs) + 1) / (nShuffle + 1);

% null hypothesis: shift is smaller in real than real
obs     = Halves_shuf.Field.real.XCor_shiftDist_sp;
shuf    = Halves_shuf.Field.shuf.HalvesXCorshiftDist_shuf;
Halves_shuf.Field.shuf.HalvesXCorshiftDist_shuf_pval    = (sum(shuf <= obs) + 1) / (nShuffle + 1);


%% shuffle for speed moduation stability/reliability

nFrames = size(smoothedSpeed, 1);
SpikeLia = false(nFrames,1);
spikeFrames_temp = spikeFrames;
spikeFrames_temp(isnan(spikeFrames_temp)) = [];
SpikeLia(spikeFrames_temp) = true;

minBinTime = 30;

%preallocate outputs
r_Pearson_shuf = nan(nShuffle,1);
r_Spearman_shuf = nan(nShuffle,1);

% ---- precompute halves indices ----
if strcmp(doCompare, 'Stability')
    idxHalve{1} = 1:midFrame-1;
    idxHalve{2} = midFrame:nFrames;
elseif strcmp(doCompare, 'Reliability')
    idxHalve{1} = find(mod(allFrame_vid,2) == 0);
    idxHalve{2} = find(mod(allFrame_vid,2) == 1);
end

parfor i_shuffle = 1:nShuffle
    % for i_shuffle = 1:nShuffle
    % --- altanative of circshift, indexing---
    shift = randi(nFrames);
    spikeIdx_orig = find(SpikeLia);
    shuffIdx = mod(spikeIdx_orig + shift - 1, nFrames) + 1;

    shuffSpikeLia = false(nFrames,1);
    shuffSpikeLia(shuffIdx) = true;

    if strcmp(doCompare, 'Stability')
        spikeIdx_halves = { find(shuffSpikeLia(1:midFrame-1)), find(shuffSpikeLia(midFrame:end)) + midFrame - 1 };
    else % Reliability
        spikeIdx_halves = { find(shuffSpikeLia(mod(allFrame_vid,2)==0)), find(shuffSpikeLia(mod(allFrame_vid,2)==1)) };
    end

    ratePerBin_halves = nan(2,20);
    local_idxHalve = idxHalve;
    local_smoothedSpeed = smoothedSpeed;
    for i = 1:2
        ratePerBin_halves(i,:) = local_ratePerBin(spikeIdx_halves{i}, local_smoothedSpeed(local_idxHalve{i}), VideoSampling, minBinTime);
    end

    [r_Pearson_shuf(i_shuffle,1), ~] = corr(ratePerBin_halves(1,:)', ratePerBin_halves(2,:)', 'Rows','complete');
    [r_Spearman_shuf(i_shuffle,1), ~] = corr(ratePerBin_halves(1,:)', ratePerBin_halves(2,:)', 'Rows','complete', 'Type', 'Spearman');

end

% toc


%% p-value
obs = Halves_shuf.Speed.real.r_Pearson;
shuf = r_Pearson_shuf;
pval_Pearson     = (sum(shuf >= obs) + 1) / (nShuffle + 1);

obs = Halves_shuf.Speed.real.r_Spearman;
shuf = r_Spearman_shuf;
pval_Spearman     = (sum(shuf >= obs) + 1) / (nShuffle + 1);


Halves_shuf.Speed.shuf.pval_Pearson = pval_Pearson;
Halves_shuf.Speed.shuf.pval_Spearman = pval_Spearman;

Halves_shuf.Speed.shuf.HalvesCorr_Pearson_shuf = r_Pearson_shuf;
Halves_shuf.Speed.shuf.HalvesCorr_Spearman_shuf = r_Spearman_shuf;

end


function ratePerBin = local_ratePerBin(spikeIdx, speed, VideoSampling, minBinTime)
edges = 0:20;
nBin = length(edges)-1;
valid = ~isnan(speed);
if isempty(spikeIdx)
    ratePerBin = nan(1,nBin);
    return
end
[N, ~] = histcounts(speed(spikeIdx(spikeIdx<=numel(speed))), edges);

timePerBin = histcounts(speed(valid), edges) / VideoSampling;
ratePerBin = N ./ (timePerBin + eps);

ratePerBin(timePerBin<minBinTime) = NaN;
end



%% Compute and show fine rate map
function [rateMap_Fine_sp, PeakFiringRate_sp, xs_sp, ys_sp, factor] = ...
    func_ComputeAndShowFineRateMap(AnimalTrack, Edges_sp, i_Cell, spikeMap_Rawcount_sp, rateMap_thr_sp, rateMapUniform_sp, binSize_cm, rateMap_Gaussian_sp)
%%
thr_sp = 3;

% ---- grid setup ----
[xDim, yDim, xPlot, yPlot, X_new, Y_new] = ...
    func_setupGrid(AnimalTrack, Edges_sp{thr_sp, i_Cell}, binSize_cm);

% ---- compute ----
rateMap_Gaussian_cell = rateMap_Gaussian_sp(:, i_Cell);
[rateMap_Fine_sp, PeakFiringRate_sp, xs_sp, ys_sp, factor] = ...
    func_computeFineRateMap(rateMap_Gaussian_cell, binSize_cm, xDim, yDim);

% ---- plot ----
func_plotFineRateMaps(spikeMap_Rawcount_sp{thr_sp,i_Cell}, rateMap_thr_sp{thr_sp,i_Cell}, rateMapUniform_sp{thr_sp,i_Cell}, ...
    rateMap_Fine_sp{thr_sp}, xPlot, yPlot, X_new, Y_new, xs_sp{thr_sp}, ys_sp{thr_sp}, xDim, yDim);

end

% helper functions
function [xDim, yDim, xPlot, yPlot, X_new, Y_new] = ...
    func_setupGrid(AnimalTrack, Edges, binSize_cm)
cropRect = AnimalTrack.TrackData.cropRect;
psh = AnimalTrack.TrackData.psh; psw = AnimalTrack.TrackData.psw;
xDim = cropRect(3)*psw; yDim = cropRect(4)*psh;

xEdges = Edges{1}; yEdges = Edges{2};
xPlot = xEdges(1:end-1); yPlot = yEdges(1:end-1);
xGrid = 0:binSize_cm:ceil(xDim/binSize_cm)*binSize_cm;
yGrid = 0:binSize_cm:ceil(yDim/binSize_cm)*binSize_cm;
[X_new, Y_new] = meshgrid(xGrid, yGrid);
end

%
function [rateMap_Fine_sp, PeakFiringRate_sp, xs_sp, ys_sp, factor] = ...
    func_computeFineRateMap(rateMap_Gaussian_cell, binSize_cm, xDim, yDim)
%%
factor = binSize_cm;
nThr = size(rateMap_Gaussian_cell,1);
rateMap_Fine_sp = cell(nThr,1);
PeakFiringRate_sp = cell(nThr,1);
xs_sp = cell(nThr,1); ys_sp = cell(nThr,1);

for k = 1:nThr
    map = rateMap_Gaussian_cell{k};
    if isempty(map), continue; end
    mapSmooth = imresize(map, factor, 'bilinear');
    xs = linspace(0, xDim, size(mapSmooth,2));
    ys = linspace(0, yDim, size(mapSmooth,1));
    rateMap_Fine_sp{k} = mapSmooth;
    PeakFiringRate_sp{k} = max(mapSmooth, [], 'all', 'omitnan');
    xs_sp{k} = xs; ys_sp{k} = ys;
end

end

%
function func_plotFineRateMaps(spkMap, rateMap_thr, rateMapUniform, ...
    mapFine, xPlot, yPlot, X_new, Y_new, xs, ys, xDim, yDim)
% colormap jet
subplot(563); func_plotMap_Fine(xPlot, yPlot, spkMap, 'Adap. bin spike map', xDim, yDim)
subplot(564); func_plotMap_Fine(xPlot, yPlot, rateMap_thr, 'Adap. bin rate map', xDim, yDim)
subplot(565); func_plotMap_Fine(X_new, Y_new, rateMapUniform, 'Rescaled rate map', xDim, yDim)
subplot(566); func_plotMap_Fine(xs, ys, mapFine, 'Smoothed fine rate map', xDim, yDim)
end

%
function func_plotMap_Fine(X, Y, Z, titleStr, xDim, yDim)
Z(isnan(Z)) = 0;
pcolor(X, Y, Z); shading flat; axis tight equal off
set(gca, 'YDir', 'reverse'); xlim([0 xDim]); ylim([0 yDim]);
title(titleStr);
end


