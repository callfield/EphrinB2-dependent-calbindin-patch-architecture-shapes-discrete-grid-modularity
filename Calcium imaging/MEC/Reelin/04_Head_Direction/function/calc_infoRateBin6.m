% Information rate using 6-degree bins
function InfoRateBin = calc_infoRateBin6(DATA, caFr)
% DATA: [angle_deg, spike_or_rate] per frame
% caFr: calcium imaging frame rate (Hz)
% Angles are binned every 6 degrees; treat -180 degrees as equivalent to 180 degrees for the wrap.

% Define 6-degree angle bin edges (from 180 down to just above -180)
angleBins = (180:-6:-179)';  % 180, 174, ..., (wrap handled below)

% Precompute totals for information-rate calculation
TotalNum = length(DATA(:,1));                     % total samples
TotalFR  = mean(DATA(:,2)) * caFr;                % overall firing rate (Hz)
tmpInfoRate = zeros(length(angleBins), 1);        % per-bin contributions

% Compute averages in each 6-degree bin.
% InfoRate contribution per bin: P(theta) * FR(theta) * log2(FR(theta)/TotalFR),
% where P(theta) is the occupancy proportion of that bin.
% Note: The final bin explicitly covers (-180, last_edge].
for i = 1:length(angleBins)
    if i < length(angleBins)
        % Bin: (angleBins(i+1), angleBins(i)]  (6 degrees wide)
        inBin = (DATA(:,1) <= angleBins(i)) & (DATA(:,1) > angleBins(i+1));
    else
        % Last bin: (-180, angleBins(i)]
        inBin = (DATA(:,1) <= angleBins(i)) & (DATA(:,1) > -180);
    end

    tmpNum = sum(inBin);
    P      = tmpNum / TotalNum;
    tmpFR  = mean(DATA(inBin, 2)) * caFr;         % mean FR in bin (Hz)

    tmpInfoRate(i,1) = P * tmpFR * log2(tmpFR / TotalFR);
end

% Sum of information-rate contributions (ignore NaNs)
InfoRateBin = sum(tmpInfoRate, 'omitnan');
end
