function [GSratio, GSpeak] = Fn_2D_gSclOri(CDir, data, samplename, s, t)

SObias = 2; % Scale : Orientation = 2:1
GSratio = nan;

%%
Gscale = round(data(:,3), 1);
Gori = round(data(:,4), 1);

gSO_map = zeros(910, 1000 * SObias); % -45 to 45 degrees, 0 to 100 cm

for c = 1:size(Gscale,1)
    tmp_gSO_map = zeros(910, 1000 * SObias); % -45 to 45 degrees, 0 to 100 cm
    Y = Gori(c,1) * 10 + 450;
    X = Gscale(c,1) * 10 * SObias;
    tmp_gSO_map(Y,X) = data(c,12);
    tmp_gSO_map = imgaussfilt(tmp_gSO_map, 10 * SObias, 'FilterSize', 30 * SObias + 1);
    tmp_gSO_map = imgaussfilt(tmp_gSO_map, 20 * SObias, 'FilterSize', 80 * SObias + 1);

    gSO_map = gSO_map + tmp_gSO_map;
end

[cent, varargout] = FastPeakFind(gSO_map, 0); % Find peaks in the scale-orientation map.

GSdist = sum(gSO_map);

clf

c = gray(255);
c(:,3) = c(:,3) * 0.1 + 0.9; c(:,1:2) = c(:,1:2) * 0.7 + 0.3;
img_gSO_map = -1 * gSO_map / max(gSO_map, [], 'all');
colormap(c)

imagesc(img_gSO_map)
hold on
scatter(Gscale * 10 * SObias, Gori * 10 + 450, 30, 'k', 'LineWidth', 2)

tmp = [];
for c = 1:length(cent) / 2
    tmp(c) = img_gSO_map(cent(2 * c), cent(2 * c - 1)) * -1;
end

GSpeak = [];
if size(tmp,2) > 1
    tmp2 = maxk(tmp, 2);
    t1 = find(tmp == tmp2(1));
    t2 = find(tmp == tmp2(2));
    plot(cent(2 * t1 - 1), cent(2 * t1), 'r*', 'MarkerSize', 7, 'LineWidth', 1.5)
    plot(cent(2 * t2 - 1), cent(2 * t2), 'r*', 'MarkerSize', 7, 'LineWidth', 1.5)
    text(cent(2 * t1 - 1), 140, num2str(round(cent(2 * t1 - 1) / (10 * SObias), 1)), 'Color', 'k', 'FontSize', 20)
    text(cent(2 * t2 - 1), 190, num2str(round(cent(2 * t2 - 1) / (10 * SObias), 1)), 'Color', 'k', 'FontSize', 20)
    xline(cent(2 * t1 - 1), '--k', 'LineWidth', 2)
    xline(cent(2 * t2 - 1), '--k', 'LineWidth', 2)

    tmp_scale = round(cent(2 * t1 - 1) / (10 * SObias), 1);
    tmp_ori = round((cent(2 * t1) - 450) / 10, 1);
    tmp_scale2 = round(cent(2 * t2 - 1) / (10 * SObias), 1);
    tmp_ori2 = round((cent(2 * t2) - 450) / 10, 1);

    GSratio = max(tmp_scale / tmp_scale2, tmp_scale2 / tmp_scale);
    text(200 * SObias, 850, strcat("GS Ratio:", num2str(round(GSratio,2))), 'Color', 'red', 'FontSize', 20)

    GSpeak(1,:) = [tmp_scale, tmp_ori];
    GSpeak(2,:) = [tmp_scale2, tmp_ori2];
    [B, I] = sort(GSpeak(:,1));
    GSpeak = GSpeak(I,:);
end

xlim([200 800] * SObias);
xticks([200 300 400 500 600 700 800] * SObias);
xticklabels([20 30 40 50 60 70 80])
ylim([100 910]);
yticks([250 460 660 860]);
yticklabels([-20 0 20 40]);
set(gca, 'YDir', 'normal')
ax = gca;
ax.XAxis.FontSize = 13
ax.YAxis.FontSize = 13
ylabel("Grid Orientation", 'Fontsize', 25)
xlabel("Grid Scale (cm)", 'Fontsize', 25)
title(strcat(samplename, " T", num2str(t)), 'FontSize', 30)

exportgraphics(gcf, strcat(CDir, "\GSOri_2D\Fig 2D ", samplename, " T", num2str(t), ".jpg"), 'Resolution', 300)
clf

end
