function [fig_field, FieldAnalysis, ACGanalysis] = func_FieldMetrics_260204_PeakRate30(binSize_cm, factor, i_Cell, rateMap_Fine_sp, thr_sp, mice_ind, Day, occMapUniform_Gaussian_sp, AnimalTrack)

%% analysis start
% mice_str = {'p13', 'p14','p15', 'p16'};
% mice_str_Folder = {'P13', 'P14','P15', 'P16'};
mice_str = {'n12', 'n14'};
mice_str_Folder = {'N12', 'N14'};


%% load data
Dir = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\Control\redefine arena border\RedefinedROIs';

fileList = dir(fullfile(Dir, '*_OFDay*.mat')); 

dayStr = num2str(Day);
pattern = ['_OF[Dd]ay', dayStr];
matchIdx = find(~cellfun(@isempty, regexp({fileList.name}, pattern)));

ArenaNewROI = load(fullfile(fileList(matchIdx(1)).folder, fileList(matchIdx(1)).name));

%% focus just on images at >2cm/s running
thrOCC = occMapUniform_Gaussian_sp{thr_sp};
BorderSize =  numel(thrOCC) * binSize_cm^2;
ActiveArea = nnz(~isnan(thrOCC)) * binSize_cm^2;

rMap_fine = rateMap_Fine_sp{thr_sp};


%% field detection
[~, ~, FieldStats] = func_FieldAnalysis_250926(rMap_fine, thr_sp, binSize_cm, factor);

% fig_field = figure;
fig_field = figure('Visible','off');

set(gcf,'Position',[-1856         494         838         353]);
tiledlayout(2,5, 'TileSpacing','compact','Padding','compact');

s_animal = sprintf( strcat('Mouse#', mice_str_Folder{mice_ind}, " Day%d"), Day);
s_CellID = sprintf('Cell ID #%d', i_Cell);
s_figtitle = strcat(s_animal, ", ", s_CellID);
sgtitle(strcat(s_figtitle, " field analysis"));

for iThr = 1:height(FieldStats)
    nexttile(iThr);
    rateMap_thr = FieldStats.thresholded_image{iThr};
    imagesc(rateMap_thr); axis tight equal off; colormap jet
    thr = FieldStats.FDthreshold(iThr);
    title(sprintf('thr: %.1f of peak', thr));
end

fprintf(strcat(s_figtitle, " field analysis\n"))

%% cut <5pix neck with waterthreshold
iThr = 3;
img = FieldStats.thresholded_image{3};

bw = img > 0;

cutpx = 5;
se = strel('disk',cutpx);
bw2 = bw;
bw2 = imerode(bw2, se);
bw2 = imdilate(bw2, se);
bw_reconstruct = imreconstruct(bw2, bw); % remove areas if imerode completely remove them

D2 = -bw2;
Ld2 = watershed(D2);
bw_cut = bw_reconstruct;
bw_cut(Ld2 == 0) = 0;

CC = bwconncomp(bw_cut, 8);
L = labelmatrix(CC); 

img_separate = bw_cut.*img;
nexttile(6);
imagesc(img_separate)
axis tight equal off; colormap jet
thr = FieldStats.FDthreshold(iThr);
title(sprintf('watershed after erode 5pix, thr: %.1f', thr));


%% peak detection
thr_val = FieldStats.threshold_val(iThr);

A1 = img_separate;
rateMapPad = padarray(A1, [1 1], 0); % zero padding

A2 = rateMapPad;
[X,Y] = meshgrid(1:size(A2,2), 1:size(A2,1));
TF2 = islocalmax2(A2, MinProminence = thr_val);

hold on
p = plot(X(TF2), Y(TF2), "kx", 'LineWidth',2, MarkerSize=8); %#ok<NASGU>
XY_peaks = [X(TF2)-1, Y(TF2)-1, A2(TF2)]; % -1: remove padding

%% draw re-defined arena ROI
NewArenaROI_raw = ArenaNewROI.RoiPos_enclose;
cropRect = ArenaNewROI.cropRect;

psh = AnimalTrack.TrackData.psh; psw = AnimalTrack.TrackData.psw;

NewArenaROI = [NewArenaROI_raw(:,1) - cropRect(1), NewArenaROI_raw(:,2) - cropRect(2)];
NewArenaROI = [NewArenaROI(:,1)*psw, NewArenaROI(:,2)*psh];
hold on
plot(NewArenaROI(:,1), NewArenaROI(:,2), 'w-', 'LineWidth',2)


%% re-analysis of fields
Length_per_bin = binSize_cm/factor; % cm
Area_per_bin   = Length_per_bin^2;

stats_re = regionprops(L, img_separate, ...
    "Area","Centroid","Circularity","MajorAxisLength","MinorAxisLength", ...
    "Orientation","MaxIntensity","MeanIntensity","WeightedCentroid","PixelIdxList");

for i = 1:numel(stats_re)
    [~, idxLocal] = max(img_separate(stats_re(i).PixelIdxList));
    maxIdx = stats_re(i).PixelIdxList(idxLocal);
    [row,col] = ind2sub(size(img_separate), maxIdx);
    stats_re(i).MaxIPos = [col,row];
    stats_re(i).Area_real = stats_re(i).Area * Area_per_bin;
end

L_temp = im2double(L)*255;
points = XY_peaks(:,[1 2]);
keepPoints = nan(size(points,1), 4);  % [x y label]
for i = 1:size(points,1)
    x = points(i,1);
    y = points(i,2);
    I = XY_peaks(i,3);
    if x >= 1 && x <= size(L,2) && y >= 1 && y <= size(L,1)
        labelVal = L_temp(y,x);
        if labelVal > 0
            % keepPoints = [keepPoints; x y labelVal];
            keepPoints(i,:) = [x, y, labelVal, I];
        end
    end
end

for i = 1:max(L(:))
    ind = find(keepPoints(:,3) == i);
    stats_re(i).LocalPeakXY = keepPoints(ind, [1, 2]);
    stats_re(i).LocalPeakI = keepPoints(ind, 4);
end

[MaxBorderScore, borderScores, allFieldsScore, Dist] = func_ComputeBorderScore(img_separate, 0);

roi = NewArenaROI;
% wallThresh = 1;
wallThresh = 5;
thr_border = 0;
[MaxBorderScore_Arena, borderScores_Arena, allFieldsScore_Arena, Dist_Arena] = func_ComputeBorderScorePolygon(img_separate, roi, thr_border, wallThresh);

for i = 1:max(L(:))
    stats_re(i).BorderScore_Old = borderScores(i);
    stats_re(i).BorderScore_Old_Dist = Dist(i);
    stats_re(i).BorderScore_PolygonArena = borderScores_Arena(i);
    stats_re(i).BorderScore_PolygonArena_Dist = Dist_Arena(i);
end

center = vertcat(stats_re.Centroid);
for i = 1:max(L(:))
    p = plot(center(i,1), center(i,2), "+", 'Color', [1 1 1]*0.7, 'LineWidth',2, MarkerSize=8 ); %#ok<NASGU>
end

%% summarize field analysis results
[PeakFR, FieldSize_Peak, FieldSize_All, FieldNum] = func_summarizeStats_Field(stats_re);

inFR = FieldStats.infield_meanfiring(iThr);
outFR = FieldStats.outfield_meanfiring(iThr);

ArenaROIs = cell2table({cropRect, NewArenaROI_raw, NewArenaROI, psw, psh}, ...
    "VariableNames",["OriginalBorder", "NewArenaROI_raw", "NewArenaROI", "psw", "psh"]);
ArenaArea = polyarea(NewArenaROI(:,1), NewArenaROI(:,2));

FieldAnalysis = {i_Cell, thr, thr_val, rMap_fine, img, img_separate, stats_re, inFR, outFR, ...
    PeakFR, FieldSize_Peak, FieldSize_All, FieldNum, ...
    MaxBorderScore, MaxBorderScore_Arena, allFieldsScore, allFieldsScore_Arena, ...
    roi, ArenaROIs, ArenaArea, BorderSize, ActiveArea,...
    FieldStats};

FieldAnalysis = cell2table(FieldAnalysis, ...
    "VariableNames",["CellID", "threshold", "threshold_val", "original_RateMap", "thres_RateMap", "watershed_RateMap", "stats","infield_meanfiring","outfield_meanfiring",...
    "PeakFR","FieldSize_Peak","FieldSize_All","FieldNum", ...
    "MaxBorderScore", "MaxBorderScore_Arena", "allFieldsScore", "allFieldsScore_Arena", ...
    "ArenaROIxy", 'ArenaROIinfo', 'ArenaArea', 'BorderSize', 'ActiveArea',...
    "AllThrImages"]);

%% show some stats
nexttile(7);

% cla
clear('s')
s{1} = sprintf('Peak firing ratio = %.3f Hz', PeakFR);
s{2} = sprintf('Peak FieldSize %.1f cm2', FieldSize_Peak);
s{3} = sprintf('FieldSizeAll = %.0f cm2', FieldSize_All);
s{4} = sprintf('FieldNum = %d', FieldNum);
s{5} = sprintf('MaxBorderScore = %.2f', MaxBorderScore_Arena);

tex = text(.0,.5, s); %#ok<NASGU>
axis off


%% Autocorrelogram analysis
%% autocorrelation with correction for zero-padding

% non-thresholded rate map
RateMap = FieldAnalysis.thres_RateMap{1};

min_overlap = 1;
Autoc = func_autocorr2_nonzeropad_250930(RateMap, min_overlap); %ACG = normxcorr2(R0, R0) make artifact due to zero-padding

% figure
nexttile(8)
imagesc(Autoc)
colormap('jet'); axis image;
clim([-1 1])
axis off
title('autocorrelogram')


%% the angle between the nearest neighbor peak (NN), the central peak (CP), and a third peak (TP)
% https://www.jneurosci.org/content/28/44/11250#sec-2

[allPeaks, CP, NN, TP, resultsTable, distCPNN, angleDeg, distCPTP] = func_ACGanalysis(RateMap, Autoc);

% Visualization
nexttile(9)
imagesc(Autoc);
colormap('jet'); axis image; %colorbar;
clim([-1 1])
hold on;

% plot peaks
plot(allPeaks(:,1), allPeaks(:,2), 'ko', 'MarkerSize',8, 'LineWidth',2);

if ~isnan(NN)
    plot([CP(1), NN(1)], [CP(2), NN(2)], 'k-', 'LineWidth',2);
end
if ~isnan(TP)
    % highlight main axis
    plot([CP(1), TP(1)], [CP(2), TP(2)], 'k-', 'LineWidth',2);
end

title(sprintf('autocorrelogram with peaks'));
xlabel('X'); ylabel('Y');
hold off;
box off
axis off

nexttile(10);
cla
clear('s')
s{1} = sprintf('distCPNN = %.1f cm', distCPNN);
s{2} = sprintf('angleDeg %.1f °', angleDeg);

tex = text(.0,.5, s); %#ok<NASGU>
axis off

%% save autocorrelogram analysis
ACGanalysis.Autoc = Autoc;
ACGanalysis.CP = CP;
ACGanalysis.NN = NN;
ACGanalysis.TP = TP;
ACGanalysis.TP_R = resultsTable;

ACGanalysis.allPeaks = allPeaks;
ACGanalysis.distCPNN = distCPNN;
ACGanalysis.distCPTP = distCPTP;
ACGanalysis.angleDeg = angleDeg;

end


%% functions
%% Place field size and other parameters
function [fig_field, FieldStats_ratio, FieldStats] = func_FieldAnalysis_250926(rateMap_Fine_sp, ~, binSize_cm, factor)
%%
rateSmooth = rateMap_Fine_sp;
Length_per_bin = binSize_cm/factor; % cm
Area_per_bin   = Length_per_bin^2;

rateSmooth_NaN = rateSmooth;
rateSmooth_NaN(rateSmooth_NaN == 0) = NaN;

fig_field = [];
FieldStats_ratio = [];

thr_field_list = [0, 0.20, 0.3, 0.5, 0.7];
FieldStats = func_runFieldAnalysis(rateSmooth, rateSmooth_NaN, ...
    thr_field_list, Area_per_bin, "ratio");
end

% --- helper functions ---
function FieldStats = func_runFieldAnalysis(rateSmooth, rateSmooth_NaN, thr_list, Area_per_bin, mode)
nThr = numel(thr_list);
FieldStats = cell(nThr, 13);

for iThr = 1:nThr
    % --- threshold
    if mode == "ratio"
        peakRate = max(rateSmooth(:), [], 'omitnan');
        thr_val  = thr_list(iThr) * peakRate;
    elseif mode == "mad" % MAD
        thr_val = std(rateSmooth(:)) * thr_list(iThr) + mean(rateSmooth(:));
    end

    [rateMap_thr, stats, inFR, outFR] = func_thresholdAndExtract(rateSmooth, rateSmooth_NaN, thr_val, Area_per_bin);
    [~, fieldScores, allFieldsScore, Dist] = func_ComputeBorderScore(rateMap_thr, thr_val);
    [PeakFR, FieldSize_Peak, FieldSize_All, FieldNum] = func_summarizeStats_Field(stats);
    FieldStats(iThr,:) = {thr_list(iThr), thr_val, rateMap_thr, stats, inFR, outFR, ...
        PeakFR, FieldSize_Peak, FieldSize_All, FieldNum, fieldScores, Dist, allFieldsScore};
end

FieldStats = cell2table(FieldStats, ...
    "VariableNames",["FDthreshold", "threshold_val", "thresholded_image","stats","infield_meanfiring","outfield_meanfiring",...
    "PeakFR","FieldSize_Peak","FieldSize_All","FieldNum", "fieldScores", "DistanceFromBorder_px", "allFieldsScore"]);
end

function [PeakFR, FieldSize_Peak, FieldSize_All, FieldNum] = func_summarizeStats_Field(stats)
if ~isempty(stats)
    [PeakFR,I]   = max([stats.MaxIntensity]);
    FieldSize_Peak = stats(I).Area_real;
    FieldNum     = numel(stats);
    FieldSize_All = sum([stats.Area_real]);
else
    PeakFR = NaN; FieldSize_Peak = NaN; FieldNum = NaN; FieldSize_All = NaN;
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
fieldMask = rateMap > thr;
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

%%

function [borderScore, fieldScores, allFieldsScore, Dist] = func_ComputeBorderScorePolygon(img_separate, roi, thr_border, wallThresh)
% ComputeBorderScorePolygon - border score from 2D rate map with arbitrary polygon ROI
%
% INPUTS:
%   rateMap    : 2D matrix (NaN for unvisited bins, numeric for firing rates)
%   roi        : ROI object from drawpolygon (with roi.Position)
%   thr        : threshold for field detection (default = 0.2 * peak)
%   wallThresh : distance [pixels] to consider "touching" the wall (default = 1)
%
% OUTPUTS:
%   borderScore : max border score among detected fields
%   fieldScores : border score for each individual field
%%

rateMap = img_separate;
if nargin < 3 || isempty(thr_border)
    thr_border = 0.2 * max(rateMap(:)); % 20% peak threshold
end
if nargin < 4 || isempty(wallThresh)
    wallThresh = 1; % within 1 pixel of wall counts as touching
end

%%

% --- polygon mask from ROI ---
polyMask = poly2mask(roi(:,1), roi(:,2), size(rateMap,1), size(rateMap,2)); %#ok<NASGU>

% --- thresholding to detect fields ---
% fieldMask = (rateMap > thr_border) & polyMask;
fieldMask = (rateMap > thr_border);
CC = bwconncomp(fieldMask);
stats = regionprops(CC, 'PixelIdxList', 'Centroid');

if isempty(stats)
    borderScore = NaN;
    fieldScores = NaN;
    allFieldsScore = NaN;
    Dist = NaN;
    return;
end

% polygon vertices
polyX = roi(:,1);
polyY = roi(:,2);

[nRows, nCols] = size(rateMap);
fieldScores = nan(1, numel(stats));
Dist = nan(1, numel(stats));
for k = 1:numel(stats)
    idx = stats(k).PixelIdxList;
    [r, c] = ind2sub(size(rateMap), idx);

    % --- centroid to boundary distance (dc) ---
    centroid = stats(k).Centroid; % (x,y)
    dc = point2polyDist(centroid(1), centroid(2), polyX, polyY);

    % --- fraction of field pixels touching wall (df) ---
    dPix = arrayfun(@(ii) point2polyDist(c(ii), r(ii), polyX, polyY), 1:numel(r));
    df = sum(dPix <= wallThresh) / numel(r);

    % --- border score ---
    fieldScores(k) = (df - (dc / max(size(rateMap)))) / (df + (dc / max(size(rateMap))));

    Dist(k) = dc;

    % disp(dc)
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


function d = point2polyDist(x, y, polyX, polyY)
% point2polyDist - shortest distance from (x,y) to polygon boundary
[~, d] = distance2curve([polyX polyY], [x y], 'linear');
end

function [pout,dist,t] = distance2curve(curve,p,method)
% DISTANCE2CURVE find closest point on curve to a given point
%
% [POUT,DIST,T] = DISTANCE2CURVE(CURVE,P,METHOD)
%
% INPUT:
%   CURVE: Nx2 array of curve vertices [x,y]
%   P    : point coordinates [x,y]
%   METHOD: 'linear' (default) or 'pchip' interpolation of curve
%
% OUTPUT:
%   POUT : closest point [x,y] on curve
%   DIST : distance from P to POUT
%   T    : parametric location along curve (0 to 1)
%
% Example:
%   t = linspace(0,2*pi,100);
%   curve = [cos(t(:)) sin(t(:))];
%   p = [0.5 0.2];
%   [pout,dist] = distance2curve(curve,p);
%   plot(curve(:,1),curve(:,2),'b-'); hold on;
%   plot(p(1),p(2),'ro',pout(1),pout(2),'kx-');
%
% Reference:
%   (C) 2011 Sven Holcombe, File Exchange #34869

if nargin<3
    method = 'linear';
end

% parametric curve length
n = size(curve,1);
tvec = linspace(0,1,n);

% interpolate curve
xi = @(tt) interp1(tvec,curve(:,1),tt,method);
yi = @(tt) interp1(tvec,curve(:,2),tt,method);

% distance function
distfun = @(tt) hypot(xi(tt)-p(1), yi(tt)-p(2));

% coarse search
tt = linspace(0,1,200);
[~,ind] = min(distfun(tt));
t0 = tt(ind);

% fminsearch refinement
t = fminsearch(distfun, t0);

% outputs
pout = [xi(t), yi(t)];
dist = distfun(t);

end


%% autocorrelogram analysis

function [allPeaks, CP, NN, TP, resultsTable, distCPNN, angleDeg, distCPTP] = func_ACGanalysis(RateMap, Autoc)
%%

% Center coordinates
R0 = RateMap;
[rows, cols] = size(R0);
centerY = rows;
centerX = cols;


% Find all local maxima above threshold
threshold_minProminence = 0.1;
threshold_minR = 0.4;
[X,Y] = meshgrid(1:size(Autoc,2), 1:size(Autoc,1));
TF = islocalmax2(Autoc, MinProminence = threshold_minProminence) & Autoc >= threshold_minR;

peakX = X(TF);
peakY = Y(TF);
Ac = Autoc(TF);

allPeaks = [peakX, peakY, Ac];
CP = [centerX, centerY];


% Identify main axis (nearest peak)
D = pdist2([peakX,peakY], [centerX, centerY]);
[D_M, I] = mink(D(:,1), 4); %#ok<ASGLU>
if isscalar(I) || isempty(I) % = length ==1
    NN = NaN;
    TP = NaN;
    resultsTable = [];
    distCPNN = NaN;
    angleDeg = NaN;
    distCPTP = NaN;
    return
elseif length(I) == 2
    NN = [peakX(I(2)) ,peakY(I(2))];
    TP = NaN;
    resultsTable = [];
    distCPNN = norm(CP - NN);
    angleDeg = NaN;
    distCPTP = NaN;
    return
end

posNN = [peakX(I(2)) ,peakY(I(2))];


% Identify second axis (third peak)
NN = posNN; % Define CP and NN
[TP, R_value, resultsTable] = func_find_third_peak(peakX, peakY, CP, NN); %#ok<ASGLU>


% angle
v1 = [NN(1)-CP(1), NN(2)-CP(2)];
v2 = [TP(1)-CP(1), TP(2)-CP(2)];

cosTheta = dot(v1,v2) / (norm(v1)*norm(v2));
cosTheta = max(min(cosTheta,1),-1);
angleDeg = acosd(cosTheta);

% outputs
distCPNN = norm(CP - NN);
distCPTP = norm(CP - TP);
end

%% find third peak in autocorrelogram mP
function [TP, R_value, resultsTable] = func_find_third_peak(peakX, peakY, CP, NN)
% FIND_THIRD_PEAK - Find the third peak (TP) in an autocorrelogram
% given the central peak (CP) and the nearest neighbor peak (NN).
%
% Inputs:
%   peakX, peakY : coordinates of all detected peaks
%   centerX, centerY : coordinates of the central peak (CP)
%   posNN : coordinates of the nearest neighbor peak (NN)
%
% Outputs:
%   TP         : selected third peak coordinates
%   R_value    : corresponding ratio value for TP
%   resultsTable : table of [peakX, peakY, R] for all peaks


% All peak coordinates
allPeaks = [peakX, peakY];

% Distance between CP and NN
distCPNN = norm(CP - NN);

% Preallocate R values
R_values = nan(size(allPeaks,1),1);

% Initialize best candidates
bestTP_noncollinear = [];
bestR_noncollinear = -inf;
bestTP_all = [];
bestR_all = -inf;

% --- Loop through all candidate peaks ---
for i = 1:size(allPeaks,1)
    TP_candidate = allPeaks(i,:);

    % Skip if candidate is CP or NN itself
    if all(TP_candidate == CP) || all(TP_candidate == NN)
        continue;
    end

    % Compute distances
    distCPTP = norm(CP - TP_candidate);
    distNNTP = norm(NN - TP_candidate);

    % Compute ratio R = dist(CP,TP) / (dist(CP,NN)+dist(NN,TP))
    R = distCPTP / (distCPNN + distNNTP);

    % Save R value
    R_values(i) = R;

    % --- Collinearity check using triangle area ---
    % If area ~ 0, CP-NN-TP are collinear
    area = abs(det([NN-CP; TP_candidate-CP]))/2;

    if area > 1e-6   % Non-collinear candidate
        if R > bestR_noncollinear
            % Update best non-collinear TP
            bestR_noncollinear = R;
            bestTP_noncollinear = TP_candidate;
        elseif abs(R - bestR_noncollinear) < 1e-12
            % Tie: choose the one with smaller angle from CP-NN
            v1 = NN - CP;
            v2_old = bestTP_noncollinear - CP;
            v2_new = TP_candidate - CP;
            ang_old = acos(dot(v1, v2_old)/(norm(v1)*norm(v2_old)));
            ang_new = acos(dot(v1, v2_new)/(norm(v1)*norm(v2_new)));
            if ang_new < ang_old
                bestTP_noncollinear = TP_candidate;
            end
        end
    end

    % --- Track global best R (for fallback case) ---
    if R > bestR_all
        bestR_all = R;
        bestTP_all = TP_candidate;
    elseif abs(R - bestR_all) < 1e-12
        % Tie: again resolve by angle
        v1 = NN - CP;
        v2_old = bestTP_all - CP;
        v2_new = TP_candidate - CP;
        ang_old = acos(dot(v1, v2_old)/(norm(v1)*norm(v2_old)));
        ang_new = acos(dot(v1, v2_new)/(norm(v1)*norm(v2_new)));
        if ang_new < ang_old
            bestTP_all = TP_candidate;
        end
    end
end

% --- Final decision: prefer non-collinear ---
if ~isempty(bestTP_noncollinear)
    TP = bestTP_noncollinear;
    R_value = bestR_noncollinear;
else
    TP = bestTP_all;   % fallback if all collinear
    R_value = bestR_all;
end

% --- Save results as a table ---
resultsTable = table(peakX, peakY, R_values, ...
    'VariableNames', {'PeakX','PeakY','R'});

end



