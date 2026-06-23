function plotColor_GridScaleVsLocation_Mod1norm(norm_data, peaks, mod1, mod2, xlims, ylims)
    figure('Position', [100, 100, 350, 500]);
    clf
    hold on

    scatter(norm_data(:,3), -2050 - norm_data(:,1), 35, 'k', LineWidth=1)
    scatter(mod1(:,3), -2050 - mod1(:,1), 40, 'r', 'filled', LineWidth=1, MarkerEdgeColor='flat')
    scatter(mod2(:,3), -2050 - mod2(:,1), 40, 'c', 'filled', LineWidth=1, MarkerEdgeColor='flat')

    % Linear regression
    x = peaks.';
    y = 0:10:1000;

    % Remove NaNs
    valid_idx = ~isnan(x);
    x_valid = x(valid_idx);
    y_valid = y(valid_idx);

    p = polyfit(x_valid, y_valid, 1);
    x_fit = linspace(xlims(1), xlims(2), 100);
    y_fit = polyval(p, x_fit);
    plot(x_fit, y_fit, 'k-', 'LineWidth', 2)
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
    eqText = sprintf('y = %.2fx + %.1f\nx @ y=%.0f: %.2f', p(1), p(2), 0, x_at_target);

    x_pos = min(xlims) + 0.1 * range(x);
    y_pos = max(ylims) - 0.1 * range(y);
    text(x_pos, y_pos, eqText, 'FontSize', 10, 'Color', 'k', ...
         'FontWeight', 'bold');

    set(gca, 'YDir', 'reverse')
    ylim(ylims)
    xlim(xlims)
    ax = gca;
    ax.XAxis.FontSize = 10;
    ax.YAxis.FontSize = 10;
    ylabel("Distance from dorsal edge of MEC (um)", 'FontSize', 12)
    xlabel("Relative Grid Scale", 'FontSize', 12)
end
