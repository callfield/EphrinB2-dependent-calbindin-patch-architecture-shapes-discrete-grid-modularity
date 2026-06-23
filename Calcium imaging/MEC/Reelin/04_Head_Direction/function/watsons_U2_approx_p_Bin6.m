% Watson's U2 test for two-sample circular data after 6-degree binning.
function [pVal, U2]=watsons_U2_approx_p_Bin6(DATA)

angleBins = (180:-6:-179)'; % 6-degree bins from 180 down to -179
DATA_Bin = DATA;
for i = 1:length(angleBins)
    if i < length(angleBins)
        DATA_Bin(DATA(:,1)<=angleBins(i) & DATA(:,1) > angleBins(i+1), 1)=angleBins(i)+3;
    else
        DATA_Bin(DATA(:,1)<=angleBins(i) & DATA(:,1) > -180, 1)=angleBins(i)+3;
    end
end

headOrientations_fired = DATA_Bin(DATA(:,2)>0,1);
headOrientations_all = DATA_Bin(:,1);
if isempty(headOrientations_fired) || isempty(headOrientations_all)
    pVal = NaN;
    U2 = NaN;
    return;
end
[pVal, U2] = watsons_U2_approx_p(headOrientations_fired,headOrientations_all);

%{
% Optional bootstrap for Watson's U2 test.
nBootstrap = 1000;
boot_pvals = zeros(nBootstrap, 1);
boot_U2 = zeros(nBootstrap, 1);

for b = 1:nBootstrap
    resampled_fired = datasample(headOrientations_fired, length(headOrientations_fired)*0.8, 'Replace', true);
    resampled_all = datasample(headOrientations_all, length(headOrientations_fired)*0.8, 'Replace', true);
    [boot_pvals(b), boot_U2(b)] = watsons_U2_approx_p(resampled_fired,resampled_all);
end

CI_pval = max(prctile(boot_pvals, [2.5, 97.5]));
CI_U2 = min(prctile(boot_U2, [2.5, 97.5]));
% pVal = max(CI_pval);
% U2 = min(CI_U2);
%}
end
