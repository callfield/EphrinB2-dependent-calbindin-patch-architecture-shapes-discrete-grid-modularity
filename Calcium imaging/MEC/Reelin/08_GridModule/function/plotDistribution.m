function peak_loc = plotDistribution(data, bandwid, faceColor, lineColor)
    peak_loc = [];
    if isempty(data)
        return;
    end

    pd = fitdist(data, 'Kernel', 'BandWidth', bandwid);
    x = linspace(min(data) - 10, max(data) + 10, 100);
    y = pdf(pd, x);
    area(x, y, 'FaceColor', faceColor, 'FaceAlpha', 0.5);
    plot(x, y, 'Color', lineColor, 'LineWidth', 2);

    [~, max_idx] = max(y);
    peak_loc = x(max_idx);
end
