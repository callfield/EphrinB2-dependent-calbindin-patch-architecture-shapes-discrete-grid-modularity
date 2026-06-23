% information rate
function InfoRate=calc_infoRate(DATA, caFr)

uniqueAngles = unique(DATA(:,1));


TotalNum = size(DATA,1); 
TotalFR = mean(DATA(:,2)) * caFr; 

infoRatePerAngle = arrayfun(@(angle) ...
    (sum(DATA(:,1) == angle) / TotalNum) * ...
    (mean(DATA(DATA(:,1) == angle, 2)) * caFr) * ...
    log2((mean(DATA(DATA(:,1) == angle, 2)) * caFr) / TotalFR), ...
    uniqueAngles, 'UniformOutput', true);

InfoRate = sum(infoRatePerAngle, 'omitnan');
end