function plot_polar_spike4(data, c, score, Angstd, HeadTitle, ax)
% data: Nx2 matrix (col 1: angle in degrees, col 2: firing rate)
% c: cell index for title display
% score: [InfoRate, Watson p-value, Watson U2, Rayleigh p-value, Rayleigh z-score]
% Angstd: weighted angular standard deviation in degrees
% HeadTitle: title string
% ax: target polar axes handle

angles = deg2rad(data(:,1));   % degrees -> radians
spike_freq = data(:,2);        % firing rate

hold on;

% Remove background and axis decorations.
ax.Color = 'none';
ax.RColor = 'none';
ax.ThetaColor = 'none';
ax.GridColor = 'none';

% Draw Cartesian-like axes in polar coordinates.
max_freq = max(spike_freq) * 1.2;
if ~isfinite(max_freq) || max_freq <= 0
    max_freq = 1;
end
polarplot(ax, [0, pi], [max_freq, max_freq], 'Color', [.3 .3 .3], 'LineWidth', 1);
polarplot(ax, [pi/2, 3*pi/2], [max_freq, max_freq], 'Color', [.3 .3 .3], 'LineWidth', 1);

% Firing-rate trace.
dark_blue = [0, 0, 0.5];
p = polarplot(ax, [angles; angles(1)], [spike_freq; spike_freq(1)], '-', ...
    'Color', dark_blue, 'LineWidth', 2, 'DisplayName', 'Spike Frequency');

ax.RLim = [0, max_freq];
uistack(p, 'top');

% Build title strings from folder names.
[parentPath, currentFolder] = fileparts(pwd);
[~, parentFolder] = fileparts(parentPath);

nameTitle = strcat(parentFolder, " ", currentFolder, " cell#", num2str(c));
nameTitle = strrep(nameTitle, '_', ' ');

scoreTitle1 = strcat(num2str(round(max_freq, 2)), "Hz, InfoRate=", num2str(round(score(1), 3)));
scoreTitle2 = strcat("WatsonU2=", num2str(round(score(3), 2)), ", pVal=", num2str(round(score(2), 4)));
scoreTitle3 = strcat("RLpVal=", num2str(round(score(4), 4)), ", RL zScore=", num2str(round(score(5), 2)), ...
                     ", Angstd=", num2str(round(Angstd, 3)));

title({HeadTitle, nameTitle, scoreTitle1, scoreTitle2, scoreTitle3}, ...
      'FontSize', 13, 'Interpreter', 'none');

hold off;
end
