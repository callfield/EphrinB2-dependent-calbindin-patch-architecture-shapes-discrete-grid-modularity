function plot_GridOriVsLocation(data, xlims, ylims)
% Plot grid orientation against dorsal distance.

    hold on
    scatter(data(:,4), -2050 - data(:,1), 30, 'k', 'filled')

    set(gca, 'YDir', 'reverse')
    ylim(ylims)
    xlim(xlims)

    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    ylabel("Distance from dorsal edge of MEC (um)", 'FontSize', 12)
    xlabel("Grid Orientation (deg)", 'FontSize', 12)
end
