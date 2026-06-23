function [bestSets, bestScores] = select_gridcell_sets(Grid_Scale, DV_pos, grid_idx, ...
                                                       nCellsPerSet, nKeepSets, nRandomSamples)

% Grid_Scale, DV_pos : vectors for all cells
% grid_idx           : indices of grid cells, such as find(wNega{s,t}(:,3)>0)
% nCellsPerSet       : number of cells in one set, such as 10
% nKeepSets          : number of final sets to keep, such as 100
% nRandomSamples     : number of random trials, such as 20000

if nargin < 6
    nRandomSamples = 20000; % default
end

% Extract feature values for the target grid cells.
gs  = Grid_Scale(grid_idx);
dv  = DV_pos(grid_idx);
nGrid = numel(grid_idx);

if nGrid < nCellsPerSet
    error('Too few grid cells (nGrid=%d < nCellsPerSet=%d)', nGrid, nCellsPerSet);
end

% Z-score the feature values and make a two-dimensional feature vector.
feat = [zscore(gs), zscore(dv)];  % size: [nGrid, 2]

% Preallocate result arrays.
sets   = zeros(nRandomSamples, nCellsPerSet, 'uint16');  % indices into grid_idx
scores = nan(nRandomSamples, 1);

% Function used for distance scoring.
pairwise_score = @(F) mean(pdist(F, 'euclidean'));  % mean pairwise distance

for s = 1:nRandomSamples
    % Randomly sample nCellsPerSet cells from grid_idx.
    idx_local = randperm(nGrid, nCellsPerSet); % indices from 1..nGrid
    sets(s, :) = idx_local;

    % Evaluate distribution width in feature space.
    F = feat(idx_local, :);   % [nCellsPerSet x 2]
    scores(s) = pairwise_score(F);
end

% Sort by score in descending order.
[sortedScores, order] = sort(scores, 'descend');
sets_sorted = sets(order, :);

% Remove duplicate sets that may occur by chance.
[uniqueSets, ia] = unique(sort(sets_sorted, 2), 'rows', 'stable');
uniqueScores = sortedScores(ia);

% Keep the top nKeepSets sets.
nKeep = min(nKeepSets, size(uniqueSets, 1));
bestSets   = uniqueSets(1:nKeep, :);   % indices into grid_idx
bestScores = uniqueScores(1:nKeep);

end
