function ang_sd = weighted_angular_std_deg(direction_firing_matrix)
    % Input: direction_firing_matrix - [angles_deg, firing_rates]
    % angles_deg: angles in degrees
    % firing_rates: firing rate or spike frequency for each angle

    angles_deg = direction_firing_matrix(:, 1);
    firing_rates = direction_firing_matrix(:, 2);

    total_rate = sum(firing_rates);
    if ~isfinite(total_rate) || total_rate <= 0
        ang_sd = NaN;
        return;
    end

    theta_rad = deg2rad(angles_deg);

    sum_cos = sum(firing_rates .* cos(theta_rad));
    sum_sin = sum(firing_rates .* sin(theta_rad));
    R = sqrt(sum_cos^2 + sum_sin^2) / total_rate;
    R = min(max(R, 0), 1);
    ang_sd = sqrt(-2 * log(R)) * (180 / pi);
end
