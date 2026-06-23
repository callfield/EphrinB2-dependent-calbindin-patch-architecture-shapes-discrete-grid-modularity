function [RLpVal, RLz] = circ_rtest_Bin6(DATA)
% Rayleigh test of non-uniformity after binning angles into 6-degree bins.
% INPUT
%   DATA: [angle_deg, value] per sample
%         angle_deg in degrees (for example, -180..+180), value > 0 means "fired"
% OUTPUT
%   RLpVal: p-value of the Rayleigh test (null: uniformity on the circle)
%   RLz   : Rayleigh's Z statistic
%
% NOTE: Uses circ_rtest from the Circular Statistics Toolbox, which expects radians.

% Define 6-degree angle-bin edges.
angleBins = (180:-6:-179)';  % 180, 174, ..., -179
DATA_Bin = DATA;

% Replace each angle by the bin center for 6-degree bins.
for i = 1:length(angleBins)
    if i < length(angleBins)
        inBin = (DATA(:,1) <= angleBins(i)) & (DATA(:,1) > angleBins(i+1));
        DATA_Bin(inBin, 1) = angleBins(i) + 3;
    else
        inBin = (DATA(:,1) <= angleBins(i)) & (DATA(:,1) > -180);
        DATA_Bin(inBin, 1) = angleBins(i) + 3;
    end
end

% Use only samples with positive value, i.e. frames with spikes or firing.
headOrientations_fired = DATA_Bin(DATA_Bin(:,2) > 0, 1);

theta_rad = deg2rad(headOrientations_fired);

if isempty(theta_rad)
    RLpVal = NaN;
    RLz = NaN;
    return;
end

[RLpVal, RLz] = circ_rtest(theta_rad);
end
