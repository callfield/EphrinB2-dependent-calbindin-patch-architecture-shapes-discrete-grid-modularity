function PEAKS = ksdnorm200(data1)

    clf

    BIN = 200;
    bandwid = 0.1;
    num_subplots = length(-1950:-10:-2950);
    xlim_range = [0 3];

    ii = 1; PEAKS = [];
    for str = -1950:-10:-2950
        tmp_data1 = data1(data1(:,1) < str & data1(:,1) > str - BIN, :);
        if isempty(tmp_data1) == 1 || size(tmp_data1,1) < 3
            PEAKS(ii) = nan;
        else
            subplot(num_subplots * 2, 1, ii * 2 - 1:ii * 2);
            hold on;

            peak_loc1 = plotDistribution(tmp_data1(:,3), bandwid, [0.5 0.5 0.5], "black");

            if ~isempty(peak_loc1)
                xline(peak_loc1, '--', 'Color', 'black', 'LineWidth', 1);
            end
            PEAKS(ii) = peak_loc1;

            xlim(xlim_range);
            yticks([]);
            if ii ~= num_subplots
                xticks([]);
            else
                break
            end
        end
        ii = ii + 1;
    end
end
