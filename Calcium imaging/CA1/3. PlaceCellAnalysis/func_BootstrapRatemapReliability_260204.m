% do bootstrap for calculating CI of ratemap and speed reliability
function [fig_boot, HalvesCorr_Pearson_even, HalvesCorr_speed_Pearson_even, HalvesCorr_Pearson_Boot, CI_rateMap, CI_rateMap_90, R_boot_mean_rateMap, ...
    HalvesCorr_speed_Pearson_Boot, CI_speed, CI_speed_90, R_boot_mean_speed, HalvesCorr_Pearson_shuf, p_rank, Cliff_delta, HalvesCorr_speed_Pearson_shuf, p_rank_speed, Cliff_delta_speed] = ...
    func_BootstrapRatemapReliability_260204(thr_sp, SpeedResults, sigma, kernelSize, mice_ind, Day, i_Cell, AnimalTrack, spikeFrames, spikePos, binSize_cm)


%% analysis start
% mice_str = {'p13', 'p14','p15', 'p16'};
% mice_str_Folder = {'P13', 'P14','P15', 'P16'};
% mice_str = {'n12', 'n14'};
mice_str_Folder = {'N12', 'N14'};

s_animal = sprintf( strcat('Mouse#', mice_str_Folder{mice_ind}, " Day%d"), Day);
s_CellID = sprintf('Cell ID #%d', i_Cell);
s_figtitle = strcat(s_animal, ", ", s_CellID);

fprintf(strcat(s_figtitle, " bootstrap analysis\n"))

mouseTimePosSpeed = [AnimalTrack.all_trk(:,1) - AnimalTrack.all_trk(1,1), AnimalTrack.all_trk(:,2:3), AnimalTrack.meanVreshape'];

Time_aligned = mouseTimePosSpeed(:,1);
VideoSampling   = numel(Time_aligned) / Time_aligned(end);
mousePos = mouseTimePosSpeed(:,[2,3]);

occThresh_reliability = 1;

runFrames = SpeedResults.idx{thr_sp}; %Index of running frame, extracted from allFrame_vid

%% normal reliability computation
runFrames_halve = cell(2,1);
runFrames_halve{1}  = runFrames(mod(runFrames,2) == 0); %even frames
runFrames_halve{2} = runFrames(mod(runFrames,2) == 1);  %odd frames
% end

[HalvesCorr_Pearson_even, p_Pearson_even] = ...
    func_HalveCompare_forBootstrap( runFrames_halve, mousePos, spikeFrames, spikePos, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, occThresh_reliability);

%% normal reliability for speed
smoothedSpeed = mouseTimePosSpeed(:,4);
allFrame_vid = (1:length(smoothedSpeed))';

% calc speed modulation reliablity
spikeFrames_halve = cell(2,1);
smoothedSpeed_halve = cell(2,1);
for i = 1:2
    % elseif strcmp(doCompare, 'Reliability')
    if i == 1
        spikeFrames_halve{i}  = spikeFrames(mod(spikeFrames,2) == 0); %even frames
        smoothedSpeed_halve_temp = smoothedSpeed;
        smoothedSpeed_halve_temp(mod(allFrame_vid,2) == 1) = NaN;
    elseif i == 2
        spikeFrames_halve{i}  = spikeFrames(mod(spikeFrames,2) == 1); %odd frames
        smoothedSpeed_halve_temp = smoothedSpeed;
        smoothedSpeed_halve_temp(mod(allFrame_vid,2) == 0) = NaN;
    end
    smoothedSpeed_halve{i} = smoothedSpeed_halve_temp;
end

[HalvesCorr_speed_Pearson_even, p_speed_Pearson_even] = ...
    func_HalveCompare_speed_forBootstrap( spikeFrames_halve, smoothedSpeed_halve, VideoSampling);


%% bootstrap reliability computation

nShuffle = 1000;
HalvesCorr_Pearson_Boot = nan(nShuffle, 1);
HalvesCorr_speed_Pearson_Boot = nan(nShuffle,1);

HalvesCorr_Pearson_shuf = nan(nShuffle, 1);
HalvesCorr_speed_Pearson_shuf = nan(nShuffle, 1);

N = numel(allFrame_vid);

% tic
% for i_shuffle = 1:nShuffle
parfor i_shuffle = 1:nShuffle
    % --- Split-half bootstrap ---
    idx = randperm(N);
    firstHalfBlocks  = idx(1:floor(N/2));
    secondHalfBlocks = idx(floor(N/2)+1:end);

    local_runFrames = runFrames;
    % --- calc rate map reliablity ---
    runFrames_halve = cell(2,1);
    runFrames_halve{1} = local_runFrames(ismember(local_runFrames, firstHalfBlocks));
    runFrames_halve{2} = local_runFrames(ismember(local_runFrames, secondHalfBlocks));

    [HalvesCorr_Pearson, ~] = ...
        func_HalveCompare_forBootstrap( runFrames_halve, mousePos, spikeFrames, spikePos, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, occThresh_reliability);

    % --- save results ---
    HalvesCorr_Pearson_Boot(i_shuffle) = HalvesCorr_Pearson;

    % shuffle
    [Idx_Lia, ~] = ismember(runFrames, spikeFrames); %spike index during running
    nRunFrames = length(Idx_Lia);
    shift = randi(nRunFrames);
    shuffIdx_logical = circshift(Idx_Lia, shift); % spike timing is shifted within running duration
    spikeFrames_inRun = runFrames(shuffIdx_logical); %shuffle
    spikePos_inRun = mousePos(spikeFrames_inRun,:); %spike position at shuffle

    [HalvesCorr_Pearson_shuf(i_shuffle), ~] = ...
        func_HalveCompare_forBootstrap( runFrames_halve, mousePos, spikeFrames_inRun, spikePos_inRun, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, occThresh_reliability);

    % --- calc speed modulation reliablity ---
    spikeFrames_halve = cell(2,1);
    smoothedSpeed_halve = cell(2,1);
    spikeFrames_halve{1} = spikeFrames(ismember(spikeFrames, firstHalfBlocks));
    spikeFrames_halve{2} = spikeFrames(ismember(spikeFrames, secondHalfBlocks));

    smoothedSpeed_halve{1} = smoothedSpeed;
    smoothedSpeed_halve{1}(secondHalfBlocks) = NaN;
    smoothedSpeed_halve{2} = smoothedSpeed;
    smoothedSpeed_halve{2}(firstHalfBlocks) = NaN;

    [HalvesCorr_speed_Pearson, ~] = ...
        func_HalveCompare_speed_forBootstrap( spikeFrames_halve, smoothedSpeed_halve, VideoSampling);

    % --- save results ---
    HalvesCorr_speed_Pearson_Boot(i_shuffle) = HalvesCorr_speed_Pearson;

    % shuffle
    nFrames = size(smoothedSpeed, 1);
    shift = randi(nFrames);
    SpikeLia = false(nFrames,1);
    local_spikeFrames = spikeFrames;
    local_spikeFrames(isnan(local_spikeFrames)) = [];
    SpikeLia(local_spikeFrames) = true;
    SpikeLia_shuff = circshift(SpikeLia, shift); % spike timing is shifted within running duration
    spikeFrames_shuff = find(SpikeLia_shuff);

    spikeFrames_halve_shuff = cell(2,1);
    spikeFrames_halve_shuff{1} = spikeFrames_shuff(ismember(spikeFrames_shuff, firstHalfBlocks));
    spikeFrames_halve_shuff{2} = spikeFrames_shuff(ismember(spikeFrames_shuff, secondHalfBlocks));

    [HalvesCorr_speed_Pearson_shuf(i_shuffle), ~] = ...
        func_HalveCompare_speed_forBootstrap( spikeFrames_halve_shuff, smoothedSpeed_halve, VideoSampling);
end
% toc

% === Compute CI and mean ===
r_boot_nonNaN = HalvesCorr_Pearson_Boot(~isnan(HalvesCorr_Pearson_Boot));
CI_rateMap = prctile(r_boot_nonNaN, [5, 95]);
CI_rateMap_90 = prctile(r_boot_nonNaN, [10, 90]);
R_boot_mean_rateMap = mean(r_boot_nonNaN);

r_boot_nonNaN = HalvesCorr_speed_Pearson_Boot(~isnan(HalvesCorr_speed_Pearson_Boot));
CI_speed = prctile(r_boot_nonNaN, [5, 95]);
CI_speed_90 = prctile(r_boot_nonNaN, [10, 90]);
R_boot_mean_speed = mean(r_boot_nonNaN);


%% === Visualization ===
% close all

fig_boot = figure('Position',[-842   284   656   493],'Visible', 'off');
% fig_boot = figure('Position',[-842   284   656   493],'Visible', 'on');
% get(gcf,'Position')

edges = -1:0.02:1;
subplot(2,2,1)
h = histogram(HalvesCorr_Pearson_Boot, edges);
hold on;
yl = ylim;
plot([R_boot_mean_rateMap R_boot_mean_rateMap], yl, 'g-', 'LineWidth',1); % mean
plot([CI_rateMap(1), CI_rateMap(1)], yl, 'r-', 'LineWidth',1); % CI lower
plot([CI_rateMap_90(1), CI_rateMap_90(1)], yl, 'r--', 'LineWidth',1); % CI lower
if ~isnan(CI_rateMap(1))
    plot([1, 1]*HalvesCorr_Pearson_even, yl, 'g--', 'LineWidth',1); % split-half reliability
end
hold off;
xlabel('Rate map split-half r'); ylabel('Count');
title(sprintf('Rate map reliability (CI95%% = %.2f, 90%% = %.2f)', CI_rateMap(1), CI_rateMap_90(1)));
xlim([-0.5, 1])
box off
ax = gca;
ax.YLim(1) = 0;

subplot(2,2,2)
histogram(HalvesCorr_speed_Pearson_Boot, edges);
hold on;
yl = ylim;
plot([R_boot_mean_speed R_boot_mean_speed], yl, 'g-', 'LineWidth',1); % mean
plot([CI_speed(1) CI_speed(1)], yl, 'r-', 'LineWidth',1); % CI lower
plot([CI_speed_90(1) CI_speed_90(1)], yl, 'r--', 'LineWidth',1); % CI lower
if ~isnan(CI_speed(1))
    plot([1, 1]*HalvesCorr_speed_Pearson_even, yl, 'g--', 'LineWidth',1); % split-half reliability
end
hold off;
xlabel('Speed modulation split-half r'); ylabel('Count');
title(sprintf('Speed mod. reliability (CI95%% = %.2f, 90%% = %.2f)', CI_speed(1), CI_speed_90(1)));
xlim([-1, 1])
box off

sgtitle(strcat(s_figtitle, " Bootstrap"));


subplot(2,2,3)
hold on;
bootall = HalvesCorr_Pearson_Boot;
shuffall = HalvesCorr_Pearson_shuf;
bootall = bootall(~isnan(bootall));
shuffall = shuffall(~isnan(shuffall));

h = histogram(bootall, edges);
h = histogram(shuffall, edges);

ax = gca;
yl = ax.YLim;
plot([1, 1]*median(bootall, 'omitmissing'), yl, 'b-', 'LineWidth',1); % mean
plot([1, 1]*median(shuffall, 'omitmissing'), yl, 'r-', 'LineWidth',1); % mean
hold off;
ax.YLim = ax.YLim * 1.2;
xlabel('Rate map split-half r'); ylabel('Count');
title(sprintf('Rate map reliability (CI95%% = %.2f, 90%% = %.2f)', CI_rateMap(1), CI_rateMap_90(1)));
xlim([-0.5, 1])
box off
ax = gca;
ax.YLim(1) = 0;

if isempty(bootall) || isempty(shuffall)
    p_rank = NaN;
    Cliff_delta = NaN;
    func_showStatsStars(p_rank, Cliff_delta)
else
    [p_rank,h,stats] = ranksum(bootall, shuffall);
    Effect = meanEffectSize(bootall, shuffall, "Effect", "cliff" );
    Cliff_delta = Effect.Effect;
    func_showStatsStars(p_rank, Cliff_delta)
end


subplot(2,2,4)
hold on
bootall = HalvesCorr_speed_Pearson_Boot;
shuffall = HalvesCorr_speed_Pearson_shuf;
bootall = bootall(~isnan(bootall));
shuffall = shuffall(~isnan(shuffall));

h = histogram(bootall, edges);
h = histogram(shuffall, edges);

ax = gca;
yl = ax.YLim;
plot([1, 1]*median(bootall, 'omitmissing'), yl, 'b-', 'LineWidth',1); % mean
plot([1, 1]*median(shuffall, 'omitmissing'), yl, 'r-', 'LineWidth',1); % mean
hold off;
ax.YLim = ax.YLim * 1.2;
xlabel('Rate map split-half r'); ylabel('Count');
title(sprintf('Rate map reliability (CI95%% = %.2f, 90%% = %.2f)', CI_rateMap(1), CI_rateMap_90(1)));
xlim([-1, 1])
box off
ax = gca;
ax.YLim(1) = 0;


if isempty(bootall) || isempty(shuffall)
    p_rank_speed = NaN;
    Cliff_delta_speed = NaN;
    func_showStatsStars(p_rank_speed, Cliff_delta_speed)
else
    [p_rank_speed ,h,stats] = ranksum(bootall, shuffall);
    Effect = meanEffectSize(bootall, shuffall, "Effect", "cliff" );
    Cliff_delta_speed = Effect.Effect;
    func_showStatsStars(p_rank_speed, Cliff_delta_speed)
end


end




%% functions
%% speed modulation
function speedModResult = func_analyzeSpeedModulation_251007(spikeFrames, smoothedSpeed, VideoSampling, doShuffle, nShuffle, doFitlm, doFitLogistic, minBinTime)
%%
% spikeFrames   : [Ns x 1] indeces of spike frames
% smoothedSpeed : [Nframes x 1] animal velocity at each frame [cm/s]
% VideoSampling : video sampling rate [Hz] (3 0Hz)
% Option:
%   'Shuffling' : 'on' / 'off' (default 'off')
%   'nShuffle'  : number of shuffling (default 1000)
%
% % ==== processing option ====
% p = inputParser;
% addParameter(p, 'Shuffling', 'off', @(x) any(validatestring(x,{'on','off'})));
% addParameter(p, 'nShuffle', 1000, @(x) isnumeric(x) && isscalar(x) && x>=0);
% parse(p, varargin{:});
% doShuffle = strcmp(p.Results.Shuffling,'on');
% nShuffle  = p.Results.nShuffle;


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

spikeCounts = histcounts(spikeSpeed, edges);
ratePerBin  = spikeCounts ./ SpeedTimeHist;                        % spike count hist to hist of firing rate [Hz]
speedCenters = edges(1:end-1) + diff(edges)/2;

% firing rate ==0 means that data is missing
ratePerBin(ratePerBin==0) = NaN;
ratePerBin_original = ratePerBin;
ratePerBin(SpeedTimeHist<minBinTime) = NaN;

% correlation
[corr_r, corr_p] = corr(speedCenters', ratePerBin', 'Rows','complete');


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
    rateShuf = nan(nBins, nShuffle);
    I_secShuf = nan(nShuffle,1);
    I_spikeShuf = nan(nShuffle,1);
    parfor i_shuffle = 1:nShuffle
        shift = randi(nFrames);
        spkShift = mod(spikeFrames + shift - 1, nFrames) + 1;
        local_smoothedSpeed = smoothedSpeed;
        spdShift = local_smoothedSpeed(spkShift);
        sc = histcounts(spdShift, edges);
        rpb = sc ./ SpeedTimeHist;
        shuffleR(i_shuffle) = corr(speedCenters', rpb', 'Rows','complete');
        rateShuf(:,i_shuffle) = rpb;

        infoShuf = func_calcSpeedInformation(rpb, SpeedTimeHist);
        I_secShuf(i_shuffle) = infoShuf.I_sec;
        I_spikeShuf(i_shuffle) = infoShuf.I_spike;
    end
    % p_shuffle = mean(abs(shuffleR) >= abs(corr_r));
    p_shuffle = (sum(abs(shuffleR) >= abs(corr_r)) + 1) / (nShuffle + 1);
    shufMean = mean(rateShuf, 2, 'omitnan');
    shufStd  = std(rateShuf, 0, 2, 'omitnan');
    zscoredRate    = (ratePerBin(:) - shufMean) ./ shufStd;

    I_secShuf_p = mean(abs(I_secShuf) >= abs(speedModResult.speedInfo.I_sec));
    I_spikeShuf_p = mean(abs(I_spikeShuf) >= abs(speedModResult.speedInfo.I_spike));

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



end

% spatial information
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


%% compute Pearson R between halved spike data rate map
function [HalvesCorr_Pearson, p_Pearson] = ...
    func_HalveCompare_forBootstrap( runFrames_halve, mousePos, spikeFrames, spikePos, binSize_cm, VideoSampling, kernelSize, sigma, AnimalTrack, occThresh_reliability)
%% calculate rate maps

rateMap_sp_halve = cell(2,1);
for i = 1:2
    runSpkIdx_halve = ismember(spikeFrames, runFrames_halve{i});
    neuronPos_halve  = spikePos(runSpkIdx_halve, :);
    mousePos_run_halve  = mousePos(runFrames_halve{i}, :);

    [occMapMerged, spkMapMerged, rateMapMerged, xEdges, yEdges, nBin, Dur_Running, areaMap, rateMapMerged_ori, occMapMerged_ori] = ...
        func_adaptiveBinRateMap(mousePos_run_halve, neuronPos_halve, AnimalTrack, binSize_cm, VideoSampling, occThresh_reliability, sigma, kernelSize); %#ok<ASGLU>

    Edges = {xEdges, yEdges};
    [occMapNormalized, occMapGaussian, occMapUniform, occMapUniform_Gaussian, X_old, Y_old, X_new, Y_new, xPlot, yPlot, xDim, yDim] = ...
        computeGaussianOccMap(AnimalTrack, occMapMerged, Edges, areaMap, binSize_cm, sigma, kernelSize); %#ok<ASGLU>
    [rateMapUniform, rateMapUniform_Gaussian] = computeGaussianRateMap(rateMapMerged, X_old, Y_old, X_new, Y_new, sigma, kernelSize);

    rateMap_sp_halve{i} = rateMapUniform_Gaussian;
end


%% calc Pearson correlation between rate maps

A = rateMap_sp_halve{1}; B = rateMap_sp_halve{2};
A(isnan(A))=0; B(isnan(B))=0;

x = A(:);
y = B(:);
validIdx = ~isnan(x) & ~isnan(y);
if sum(validIdx) >= 3
    [HalvesCorr_Pearson, p_Pearson] = corr(x(validIdx), y(validIdx), 'Rows', 'complete');
else
    HalvesCorr_Pearson = NaN;
    p_Pearson = NaN;
end



end


%% compute Pearson R between halved spike data speed hist
function [HalvesCorr_speed_Pearson, p_speed_Pearson] = ...
    func_HalveCompare_speed_forBootstrap( spikeFrames_halve, smoothedSpeed_halve, VideoSampling)
%%
speedModResult_halve = cell(2,1);
ratePerBin_halves = nan(2,20);
for i = 1:2
    doShuffle = 0; nShuffle = 0; doFitlm = 0; doFitLogistic = 0;
    minBinTime = 60;
    spikeFrames_halve_temp = spikeFrames_halve{i};
    smoothedSpeed_halve_temp = smoothedSpeed_halve{i};
    speedModResult_halve{i} = func_analyzeSpeedModulation_251007(spikeFrames_halve_temp, smoothedSpeed_halve_temp, VideoSampling, doShuffle, nShuffle, doFitlm, doFitLogistic, minBinTime);
    ratePerBin = speedModResult_halve{i}.ratePerBin;
    ratePerBin_halves(i,1:20) = ratePerBin;
    speedCenters = speedModResult_halve{i}.speedCenters;

end

% if comparable bin number is only 2, R returns -1 or 1.
% Set minumum comparable number of bins
x = ratePerBin_halves(1,:)';
y = ratePerBin_halves(2,:)';
validIdx = ~isnan(x) & ~isnan(y);
if sum(validIdx) >= 3
    [r_Pearson, p_Pearson] = corr(x(validIdx), y(validIdx), 'Rows', 'complete');
else
    r_Pearson = NaN;
    p_Pearson = NaN;
end

doImShow = 0;
if doImShow
    % subplot(5,10,5)
    figure
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
        });
    xlim([0 15])
    box off
    hold on

end


HalvesCorr_speed_Pearson = r_Pearson;
p_speed_Pearson = p_Pearson;


end



%%
function func_showStatsStars(p_rank, Cliff_delta)

% marks for p values
if p_rank < 0.001
    sigMark = '***';
elseif p_rank < 0.01
    sigMark = '**';
elseif p_rank < 0.05
    sigMark = '*';
else
    sigMark = 'n.s.';
end

% marks for effect size（Cliff’s δ: small = 0.147, medium = 0.33, large = 0.47）
if abs(Cliff_delta) > 0.47
    effMark = '###';
elseif abs(Cliff_delta) > 0.33
    effMark = '##';
elseif abs(Cliff_delta) > 0.147
    effMark = '#';
else
    effMark = 'ignorable';
end


txt = {sprintf('p = %.3g  (%s)', p_rank, sigMark),...
    sprintf('Cliff\\delta = %.3f  (%s)',Cliff_delta, effMark)};

yl = ylim;
xpos = mean(xlim); 
ypos = yl(2) * 0.99; 
text(xpos, ypos, txt, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'Interpreter', 'tex');

title('Blue: bootstrap, Red: shuffle');
xlabel('Rate map split-half r');
ylabel('Count');
end

%%


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



%% Generate Gaussian Occ and RateMap

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
