function make_tag_montage(CDir, th_tags, file_pattern, out_name, nrows, ncols)
% CDir        : directory containing source images
% th_tags     : cell array of tags, such as {'r02p05', ...}
% file_pattern: sprintf pattern, such as 'MaxErrorDist_log_hist_GROUP1vsGROUP2_%s.jpg'
% out_name    : output jpg file name without extension
% nrows, ncols: number of rows and columns

    % Use a tall layout for 6 x 2 panels.
    fig_w = 600;
    fig_h = 3600;
    fig = figure('Units','pixels', ...
                 'Position',[100 100 fig_w fig_h], ...
                 'Color','w');

    % Reduce spacing between panels.
    tiledlayout(nrows, ncols, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    ntags      = numel(th_tags);
    max_panels = nrows * ncols;
    if ntags > max_panels
        warning('The number of tags (%d) exceeds the number of panels (%d). Showing the first %d tags only.',...
            ntags, max_panels, max_panels);
        ntags = max_panels;
    end

    for i = 1:ntags
        tag   = th_tags{i};
        fname = fullfile(CDir, sprintf(file_pattern, tag));
        if ~isfile(fname)
            warning('File not found: %s', fname);
            nexttile; axis off;
            continue;
        end

        I = imread(fname);

        ax = nexttile;
        imshow(I, 'Parent', ax);
        axis(ax, 'off'); hold(ax, 'on');

        % Add a bold tag label near the top of each image.
        [h, w, ~] = size(I);
        text(ax, w/4, h*0.08, tag, ...
            'FontSize', 13, ...
            'FontWeight','bold', ...
            'HorizontalAlignment','center', ...
            'Color','w', ...
            'BackgroundColor','k', ...
            'Margin',5);

        hold(ax, 'off');
    end

    % Add an overall title here if needed.
    % title_str = strrep(out_name,'_','\_');
    % sgtitle(tlo, title_str, 'FontSize', 18, 'FontWeight','bold');

    out_jpg = fullfile(CDir, [out_name '.jpg']);
    exportgraphics(fig, out_jpg, 'Resolution', 300);
    close(fig);
end
