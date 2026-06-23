function peak = plotDistribution_ksd(data, bandwid, faceColor)
    peak = [];
    if isempty(data)
        return;
    end

    pd = fitdist(data, 'Kernel', 'BandWidth', bandwid);
    x = linspace(min(data), max(data), 100);
    y = pdf(pd, x);

    area(x, y, 'FaceColor', faceColor, 'FaceAlpha', 0.5);
    hold on
    plot(x, y, 'Color', faceColor * 0.8, 'LineWidth', 2);

    [~, max_idx] = max(y);
    peak = x(max_idx);
end
