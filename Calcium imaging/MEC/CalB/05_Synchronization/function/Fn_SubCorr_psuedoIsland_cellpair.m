function [data1 data2 data3 data4] =Fn_SubCorr_psuedoIsland_cellpair( DATA, ISLAND, TMPNET, DATA2)
% Function to analyze and visualize pseudo-islands and compute ratio of significant correlations

% Initialize
clf;
NumCell = size(DATA, 1);
AA = nan(NumCell, NumCell, NumCell);
for i = 1:NumCell
    for j = i+1:NumCell
        AA(:, i, j) = 1;
    end
end

% Identify top pseudo-islands
NumPair = sum(TMPNET, 2, 'omitnan');
[B, I] = sort(NumPair, 'descend');
I = I(B > 5);
if length(I) >= 10
    I = I(1:10);
end
if isempty(I)
    data1 = [];
    data2 = [];
    data3 = [];
    data4 = [];
    return
end

mIsland = find(ISLAND == 1);
IntraIsland = cell(length(I), 1);
II = cell(length(I), 1);

for l = 1:length(I)
    IntraIsland{l} = nan(NumCell, NumCell, NumCell);
    tmp0 = find(TMPNET(:, I(l)) == 1);
    tmp = mIsland(tmp0);
    for q = 1:length(tmp)
        for qq = q+1:length(tmp)
            i = tmp(q);
            j = tmp(qq);
            IntraIsland{l}(:, i, j) = 1;
        end
    end
    II{l} = tmp;
end



% Identify trans-island
% Trans island: select the pseudo-island at least 300 um from the top reference cell.
[B, I] = sort(sum(TMPNET, 2, 'omitnan'), 'descend');
I = I(B > 5);
pIslandCellPos = DATA2(mIsland(I(1)), 1:2);
otherIslandCellPos = DATA2(mIsland, 1:2);
dists = pdist2(pIslandCellPos, otherIslandCellPos);
canditateCell = find(dists > 300);

psudo_transIdx = [];
tmp = I(ismember(I, canditateCell));
if ~isempty(tmp), psudo_transIdx = tmp(1); end

tmp0=find(TMPNET(:,psudo_transIdx)==1);
tmp=mIsland(tmp0);
trans_psudoIslandCells = tmp;

top_psudoIslandCells = II{1};




% All island
All_Island = nan(NumCell, NumCell, NumCell);
for q = 1:length(mIsland)
    for qq = 1:length(mIsland)
        i = mIsland(q);
        j = mIsland(qq);
        All_Island(:, i, j) = 1;
    end
end


% Make trans-island pair matrix
top_psudoIslandCells = II{1};
psudo_transIsland = nan(NumCell, NumCell, NumCell);
for q = 1:length(trans_psudoIslandCells)
    for qq = 1:length(top_psudoIslandCells)
        i = trans_psudoIslandCells(q);
        j = top_psudoIslandCells(qq);
        psudo_transIsland(:, i, j) = 1;
    end
end



% 1: ALL island -ALL island
% 2: Top 1 psudo Intra Island & paired psudo island
% 3: Average of Top 10 Intra psudo Island
% 4: trans psudo island



tmp = reshape(DATA(mIsland, :, :) .* All_Island(mIsland, :, :) * 100, [],1);
data1 = tmp(~isnan(tmp));

data2 = [];
tmp01 = reshape(DATA(II{1}, :, :) .* IntraIsland{1}(II{1}, :, :) * 100, [],1);
tmp02 = reshape(DATA(trans_psudoIslandCells, :, :) .* psudo_transIsland(trans_psudoIslandCells, :, :) * 100, [], 1);
tmp= [tmp01; tmp02];
data2 = tmp(~isnan(tmp));



data3 = [];
tmp=[];
for l = 1:min(10, numel(II))
    tmp0 = reshape(DATA(II{l}, :, :) .* IntraIsland{l}(II{l}, :, :) * 100, [],1);
    tmp = [tmp; tmp0];
end
data3 = tmp(~isnan(tmp));

tmp = reshape(DATA([top_psudoIslandCells; trans_psudoIslandCells], :, :) .* psudo_transIsland([top_psudoIslandCells; trans_psudoIslandCells], :, :) * 100, [], 1);
data4 = tmp(~isnan(tmp));



end
