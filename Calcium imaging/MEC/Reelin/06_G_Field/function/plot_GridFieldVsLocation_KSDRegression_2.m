function [slope, x_at_target] = plot_GridFieldVsLocation_KSDRegression_2(data, peaks, xlims, ylims)
% Plot grid field against dorsal distance and fit a regression line through KSD peak locations.

    hold on

    scatter(data(:,5), -2050 - data(:,1), 30, 'k', 'filled')

    x = peaks.';
    y = 100:200:900;

    valid_idx = isfinite(x);
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);

    if numel(x_valid) < 2
        error('At least two valid KSD peak locations are required for regression.');
    end

    p = polyfit(x_valid, y_valid, 1);  % Fit y = p(1)*x + p(2)
    x_fit = linspace(xlims(1), xlims(2), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'r-', 'LineWidth', 2)
    slope = p(1);

    mdl = fitlm(x_valid, y_valid);
    coeffs = mdl.Coefficients;

    p_value = coeffs.pValue(2);
    fprintf('slope p = %.4f\n', p_value)

    sx = std(x_valid);
    sy = std(y_valid);

    r_from_slope = slope * sx / sy;
    fprintf('r = %.4f\n', r_from_slope)

    x_at_target = (0 - p(2)) / p(1);
    eqText = sprintf('y = %.2fx + %.1f', p(1), p(2));

    x_pos = min(xlims) + 0.1 * range(xlims);
    y_pos = max(ylims) - 0.1 * range(ylims);
    text(x_pos, y_pos, eqText, 'FontSize', 10, 'Color', 'k', ...
         'FontWeight', 'bold');

    set(gca, 'YDir', 'reverse')
    ylim(ylims)
    xlim(xlims)
    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    ylabel("Distance from dorsal edge of MEC (um)", 'FontSize', 12)
    xlabel("Grid Field (cm)", 'FontSize', 12)
end
