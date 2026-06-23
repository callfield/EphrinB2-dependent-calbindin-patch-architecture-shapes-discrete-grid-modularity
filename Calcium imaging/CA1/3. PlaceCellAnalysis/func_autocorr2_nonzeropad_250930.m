function Autoc = func_autocorr2_nonzeropad_250930(GSrate_map, min_overlap)
    % AUTOCORR_RATE_MAP_FFT Compute 2D autocorrelogram of a rate map with overlap correction
    % GSrate_map : input rate map (2D matrix)
    % min_overlap: minimum number of overlapping pixels to compute correlation
    %
    % Returns Autoc : normalized autocorrelation map, same size as full xcorr2 output

    if nargin < 2, min_overlap = 1; end
    
    GSrate_map = double(GSrate_map); % Convert input to double for computation

    [h, w] = size(GSrate_map);       % Get height and width of the rate map

    % --- Compute numerator (cross-product) using full 2D correlation ---
    num = xcorr2(GSrate_map, GSrate_map);  % Equivalent to sum(r1.*r2) for all shifts

    % --- Compute overlap counts for each shift (to correct zero-padding) ---
    overlap = xcorr2(ones(h,w), ones(h,w)); % Number of overlapping pixels at each shift

    % --- Compute sums and squared sums for normalization ---
    sum1   = xcorr2(GSrate_map, ones(h,w));       % sum of r1 over overlapping region
    sum2   = xcorr2(ones(h,w), GSrate_map);       % sum of r2 over overlapping region
    sumsq1 = xcorr2(GSrate_map.^2, ones(h,w));    % sum of r1 squared over overlapping region
    sumsq2 = xcorr2(ones(h,w), GSrate_map.^2);    % sum of r2 squared over overlapping region

    % --- Assign number of overlapping pixels to n ---
    n = overlap;

    % --- Compute denominators for Pearson-like correlation ---
    denom1 = n .* sumsq1 - sum1.^2; % variance term for r1
    denom2 = n .* sumsq2 - sum2.^2; % variance term for r2

    % --- Compute numerator adjusted by means ---
    numer = n .* num - sum1 .* sum2; % covariance term

    % --- Initialize output and apply normalization ---
    Autoc = zeros(size(num));                     % Preallocate autocorrelation map
    valid = (n >= min_overlap) & (denom1 > 0) & (denom2 > 0); % Only valid shifts
    Autoc(valid) = numer(valid) ./ sqrt(denom1(valid) .* denom2(valid)); % Normalized correlation
end


