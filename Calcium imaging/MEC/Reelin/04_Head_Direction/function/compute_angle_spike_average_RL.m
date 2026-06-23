function [angleSpikeAvg, angleBinSpikeAvg, HD_Score_Bin1, HD_Score_Bin6] = compute_angle_spike_average_RL(DATA, caFr)
% DATA: matrix with columns [round(angle_deg), spike]
% Treat -180 as 180 for circular wrapping.
DATA(DATA(:,1) == -180, 1) = 180;

% Ensure coverage of all integer angles from -179 to 180.
allAngles = (-179:180)';

% Compute per-degree averages and convert spike probability to Hz.
uniqueAngles = unique(DATA(:,1));
spikeMeans = zeros(size(allAngles));

for i = 1:length(uniqueAngles)
    angle = uniqueAngles(i);
    spikeMeans(allAngles == angle) = mean(DATA(DATA(:,1) == angle, 2)) * caFr;
end

% Z-score across angles. If variance is zero, return zeros.
mu = mean(spikeMeans);
sigma = std(spikeMeans);
if sigma > 0
    spikeMeansZ = (spikeMeans - mu) / sigma;
else
    spikeMeansZ = zeros(size(spikeMeans));
end

% Output: [angle, mean_FR_Hz, zscore].
angleSpikeAvg = [allAngles, spikeMeans, spikeMeansZ];

% Information rate and circular tests using per-degree data.
InfoRate = calc_infoRate(DATA, caFr);
headOrientations_fired = DATA(DATA(:,2) > 0, 1);
headOrientations_all   = DATA(:,1);
if isempty(headOrientations_fired)
    pVal = NaN;
    U2 = NaN;
    RLpVal = NaN;
    RLz = NaN;
else
    [pVal, U2] = watsons_U2_approx_p(headOrientations_fired, headOrientations_all);

    % Rayleigh test (H0: uniformity; H1: unimodal directionality).
    theta_rad = deg2rad(headOrientations_fired);
    [RLpVal, RLz] = circ_rtest(theta_rad);
end

HD_Score_Bin1 = [InfoRate, pVal, U2, RLpVal, RLz];

% 6-degree binning.
angleBins = (180:-6:-179)';     % 180, 174, ..., -179
spikeMeanBins  = zeros(size(angleBins));
spikeMeanZBins = zeros(size(angleBins));

% Average spikes within each 6-degree interval.
for i = 1:length(angleBins)
    if i < length(angleBins)
        inBin = (DATA(:,1) <= angleBins(i)) & (DATA(:,1) > angleBins(i+1));
    else
        inBin = (DATA(:,1) <= angleBins(i)) & (DATA(:,1) > -180);
    end
    spikeMeanBins(i) = mean(DATA(inBin, 2)) * caFr;
end

% Z-score across 6-degree bins. If variance is zero, return zeros.
mu = mean(spikeMeanBins);
sigma = std(spikeMeanBins);
if sigma > 0
    spikeMeanZBins = (spikeMeanBins - mu) / sigma;
else
    spikeMeanZBins = zeros(size(spikeMeanBins));
end

% Output per 6-degree bin: [bin_edge_deg, mean_FR_Hz, zscore].
angleBinSpikeAvg = [angleBins, spikeMeanBins, spikeMeanZBins];

% Information rate and circular tests using 6-degree binned data.
InfoRateBin = calc_infoRateBin6(DATA, caFr);
[pVal_Bin, U2_Bin] = watsons_U2_approx_p_Bin6(DATA);

% Rayleigh test after 6-degree binning.
[RLpVal_Bin, RLz_Bin] = circ_rtest_Bin6(DATA);
HD_Score_Bin6 = [InfoRateBin, pVal_Bin, U2_Bin, RLpVal_Bin, RLz_Bin];
end
