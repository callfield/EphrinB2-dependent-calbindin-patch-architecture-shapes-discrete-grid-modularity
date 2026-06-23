function [Mod, ModScale, ModOri] = Fn_Define_2DModule(CDir, data, samplename, s, t, GSpeak)

SObias = 2; % Scale : Orientation = 2:1

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

%%
ModOri = []; ModScale = []; Mod = []; ModpeakGs = nan(2,1);
if size(GSpeak,1) > 1
    for m = 1:2
        tmp_Gscl = data(:,3);
        tmp_Gori = data(:,4);

        % Define Grid orientation of module 1.
        xx = gSO_map(:, ((GSpeak(m,1) - 5) * 10 * SObias):((GSpeak(m,1) + 5) * 10 * SObias));
        xxx = sum(xx,2);
        [pks, locs, w, p, wxPk] = findpeaks_ho(xxx, 'SortStr', 'descend', 'NPeaks', 1);
        tmp = ((round(wxPk,0) - 450) / 10 - GSpeak(m,2));
        ModOri(m,:) = sort(round((wxPk - 450) / 10 + tmp, 0));

        % Define Grid Scale of module 1.
        xx = gSO_map((ModOri(m,1) * 10 + 450):(ModOri(m,2) * 10 + 450), :);
        xxx = sum(xx);
        T = islocalmin(xxx);
        yy = find(T == 1).';
        tmp_min = max(yy(find(yy < GSpeak(m,1) * 10 * SObias)));
        tmp_max = min(yy(find(yy > GSpeak(m,1) * 10 * SObias)));
        if isempty(tmp_min)
            ModScale(m,2) = round(tmp_max / (10 * SObias), 1);
            ModScale(m,1) = GSpeak(m,1) - abs(GSpeak(m,1) - ModScale(m,2));
        else
            ModScale(m,1) = round(tmp_min / (10 * SObias), 1);
        end

        if isempty(tmp_max)
            ModScale(m,1) = round(tmp_min / (10 * SObias), 1);
            ModScale(m,2) = GSpeak(m,1) + abs(GSpeak(m,1) - ModScale(m,1));
        else
            ModScale(m,2) = round(tmp_max / (10 * SObias), 1);
        end

        Mod{m} = find(ModOri(m,1) < tmp_Gori & ModOri(m,2) > tmp_Gori & ...
            ModScale(m,1) < tmp_Gscl & ModScale(m,2) > tmp_Gscl);
    end

    % Remove duplicate cells.
    commonIdx = intersect(Mod{1}, Mod{2});
    Mod{1} = setdiff(Mod{1}, commonIdx);
    Mod{2} = setdiff(Mod{2}, commonIdx);

    ModpeakGs = nan(2,1);
    ModpeakGs(1) = gSO_map(GSpeak(1,2) * 10 + 450, GSpeak(1,1) * 10 * SObias);
    ModpeakGs(2) = gSO_map(GSpeak(2,2) * 10 + 450, GSpeak(2,1) * 10 * SObias);

    clf
    subplot(4,4,2)
    scatter(data(:,3), data(:,4), 25, 'k')
    hold on
    scatter(data(Mod{1},3), data(Mod{1},4), 20, 'r', 'filled')
    scatter(data(Mod{2},3), data(Mod{2},4), 20, 'c', 'filled')
    xticks([]);
    ax = gca; ax.YAxis.FontSize = 10

    subplot(4,4,[6 10 14])
    scatter(data(:,3), data(:,1), 25, 'k')
    hold on
    scatter(data(Mod{1},3), data(Mod{1},1), 20, 'r', 'filled')
    scatter(data(Mod{2},3), data(Mod{2},1), 20, 'c', 'filled')
    set(gca, 'YDir', 'normal')
    yticks([]); ylim([-2900 -1950])
    ax = gca;
    ax.XAxis.FontSize = 10
    ax.YAxis.FontSize = 10
    ylabel("Cell location (Ventral - Dorsal)", 'Fontsize', 12)
    xlabel("Grid Scale (cm)", 'Fontsize', 12)

    subplot(4,4,[7:8 11:12 15:16])
    scatter(data(:,2), data(:,1), 25, 'k')
    hold on
    scatter(data(Mod{1},2), data(Mod{1},1), 20, 'r', 'filled')
    scatter(data(Mod{2},2), data(Mod{2},1), 20, 'c', 'filled')
    title({strcat(samplename, " T", num2str(t)); "FoV"}, 'FontSize', 12)
    set(gca, 'YDir', 'normal')
    xticks([]); yticks([]);
    xlim([-200 720]); ylim([-2900 -1950])
    ax = gca;
    ax.XAxis.FontSize = 10
    ax.YAxis.FontSize = 10
    xlabel("Cell location (Medial - Lateral)", 'Fontsize', 12)
    legend({'Grid', 'Module 1', 'Module 2'}, 'Location', 'best')

    exportgraphics(gcf, strcat(CDir, "\Gmodule_2D_220822\Gmod 2D ", samplename, " T", num2str(t), ".jpg"), 'Resolution', 300)
    clf
end

end
