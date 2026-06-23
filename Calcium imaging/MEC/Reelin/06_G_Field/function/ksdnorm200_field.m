function PEAKS = ksdnorm200_field(data1)

    clf

    binWidth = 200;
    bandwid = 3.5;
    numSubplots = 5;
    xlimRange = [10 100];

    PEAKS = nan(1, numSubplots);

    for ii = 1:numSubplots
        regionStart = -2050 - (ii - 1) * binWidth;
        tmp_data1 = data1(data1(:,1) < regionStart & data1(:,1) > regionStart - binWidth, :);

        if isempty(tmp_data1) || size(tmp_data1, 1) < 3
            continue
        end

        subplot(numSubplots * 2, 1, ii * 2 - 1:ii * 2);
        hold on

        peak_loc1 = plotDistribution(tmp_data1(:,5), bandwid, [0.5 0.5 0.5], "black");

        if ~isempty(peak_loc1)
            xline(peak_loc1, '--', 'Color', 'black', 'LineWidth', 1);
            PEAKS(ii) = peak_loc1;
        end

        xlim(xlimRange);
        yticks([]);
        if ii ~= numSubplots
            xticks([]);
        end
    end
end
