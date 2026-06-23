clear
close all
%% load exported trace data

% mice_ind = 1:10;

% mice_str = {'n5', 'n6', 'n7', 'n8', 'n9', 'n10', 'p1', 'p7','p8', 'p11'};
% mice_str_Folder = {'N5', 'N6', 'N7', 'N8', 'N9', 'N10', 'P1', 'P7','P8', 'P11'};
% mice_str = {'p13', 'p14', 'p15', 'p16'};
% mice_str_Folder = {'P13', 'P14','P15', 'P16'};
mice_str = {'n12', 'n14'};
mice_str_Folder = {'N12', 'N14'};



for mice_ind = 1 %1:2
    for Day = 1:6


        fprintf(strcat('Mouse', mice_str_Folder{mice_ind}, " Day%d\n"), Day);
        CaFolderName = strcat('CalB', mice_str_Folder{mice_ind});


        caPath = fullfile('K:\experiments\260121 Ca imaging\analysis_PCAICA', CaFolderName, 'export');
        caFile = strcat(mice_str{mice_ind}, '_d', num2str(Day), '.csv');


        if exist(fullfile(caPath, caFile), 'file') == 0
            fprintf(strcat("File not found: ", 'Mouse', mice_str_Folder{mice_ind}, " Day%d\n"), Day);
            continue;
        end

        calciumTraces = readmatrix(fullfile(caPath, caFile));  % [frames x ROIs]
        calciumTraces_table = readtable(fullfile(caPath, caFile));  % [frames x ROIs]

        [~, name, ext] = fileparts(caFile);
        caFile_prop = strcat(name, "-props", ext);
        calciumTraces_prop = readmatrix(fullfile(caPath, caFile_prop));
        calciumTraces_prop_table = readtable(fullfile(caPath, caFile_prop));



        %% distance between cells

        dist_thr = 10; % [pix] ~= 20um

        data_prop = calciumTraces_prop;
        coordinates = data_prop(:, [6 7]);

        D = pdist2(coordinates, coordinates);
        D(D>dist_thr) = 0;


        %% load data
 
        data_dir = fullfile('K:\experiments\260121 Ca imaging\analysis_PCAICA', CaFolderName, strcat(CaFolderName, '_data'));

        patterns = { ...
            strcat(CaFolderName, '_OFDay', num2str(Day)), ...
            strcat(CaFolderName, '_Day',  num2str(Day)), ...
            strcat(CaFolderName, '_day',  num2str(Day)), ...
             };

        fnamelist = [];
        for i = 1:numel(patterns)
            try
                fnamelist = filenamelisting(data_dir, patterns{i});
                if ~isempty(fnamelist)
                    str_rule = patterns{i};
                    break;  % stop once a valid file list is found
                end
            catch
                % do nothing, try next pattern
            end
        end

        if isempty(fnamelist)
            error('No matching files found for CaFolderName = %s, Day = %d', CaFolderName, Day);
        end


        % load PCA-ICA isxd file
        clear('fname')
        A = fnamelist;
        str_rule = 'PCA-ICA';
        for i = 1:length(A)
            if  (isempty( regexp(A{i}, str_rule, 'once')) == 0)
                fname{i,1} = A{i,1};
                fname{i,2} = A{i,2};
                fname{i,3} = A{i,3};
            end
        end
        x = find(cellfun('isempty',fname(:,1)));
        fname(x,:) = [];
        fnamelist_ICA = fname;

        [~, I] = max( cell2mat(fnamelist_ICA(:,3)));
        PCAICAfile = fullfile(data_dir, fnamelist_ICA{I,1});
        cell_set = isx.CellSet.read(PCAICAfile);


        % load max projection isxd file
        clear('fname')
        str_rule = 'Maximum Image';
        for i = 1:length(A)
            if  (isempty( regexp(A{i}, str_rule, 'once')) == 0)
                fname{i,1} = A{i,1};
                fname{i,2} = A{i,2};
                fname{i,3} = A{i,3};
                fname{i,3} = A{i,4};
            end
        end
        x = find(cellfun('isempty',fname(:,1)));
        fname(x,:) = [];
        fnamelist_MP = fname;
        [~, I] = max( cell2mat(fnamelist_MP(:,3)));
        mpim  = fullfile(data_dir, fnamelist_MP{I,1});

        image = isx.Image.read(mpim);
        image_data = image.get_data();


        figure
        imshow(( image_data./prctile(image_data(:), 95) ) * 0.6)
        hold on
        for kk = 1:cell_set.num_cells
            c = coordinates(kk,:);
            text(c(1), c(2)-3, sprintf('%d', kk), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'Color', 'w');
        end
        hold off


        %% show image
        % accumulate cell masks & boundaries in one loop
        nCells = cell_set.num_cells;
        CellImage = zeros(size(image_data));
        Boundaries = cell(nCells,1);

        for kk = 1:nCells
            img = cell_set.get_cell_image_data(kk - 1);
            CellMask = func_findboundary(img);
            B = bwboundaries(CellMask);

            if ~isempty(B)
                Boundaries{kk} = B{1}; % cell boundaries
            end

            % mask with cell indices
            CellImage(CellMask > 0) = kk;
        end

        % prepare base image (normalized)
        A = ( image_data ./ prctile(image_data(:),95) ) * 0.6;

        % show overlay
        figure;
        I3 = labeloverlay(A, CellImage, 'Colormap','jet','Transparency',0.4);
        imshow(I3); hold on;

        % plot all boundaries
        for kk = 1:nCells
            if ~isempty(Boundaries{kk})
                plot(Boundaries{kk}(:,2), Boundaries{kk}(:,1), 'r', 'LineWidth',0.5);
            end
        end


        %% read traces

        close all

        i = 32;
        Time = calciumTraces(3:end,1);
        tracedata = calciumTraces(3:end,2:end);
        Trace = tracedata(:,i);
        Max = max(tracedata(:));

        % figure
        % for i = 1:100 %size(tracedata,2)
        %     Trace = tracedata(:,i);
        %     plot(Time, Trace/Max - (i-1))
        %     % plot(Time, Trace/Max)
        %     hold on
        % end


        %% cell validation

        close all

        % shared variables
        Time = calciumTraces(3:end,1);
        tracedata = calciumTraces(3:end,2:end);  % [frames × cells]
        nCells = size(tracedata,2);

        Diff_time = diff(Time);
        Fs = 1 / Diff_time(1);
        full_ts = Time(1):1/Fs:Time(end);

        % interporate for concatenated data
        interp_traces = interp1(Time, tracedata, full_ts, 'linear'); % interp for each colomn, [length(full_ts) × nCells]

        % output variables
        Metric = zeros(nCells,4);
        PeakIdx = cell(nCells,2);
        trace_smooth_all = zeros(nCells, length(full_ts));

        % compute metrics
        window_size = 500; % baseline drift correction
        for i = 1:nCells %470 %4 %1:nCells

            interp_trace = interp_traces(:,i);

            % baseline drift correction
            baseline_drift = movmedian(interp_trace, window_size);
            trace_detrended = interp_trace - baseline_drift;

            % smoothing
            trace_smooth = smoothdata(trace_detrended, 'movmean', 5);
            trace_smooth_all(i,:) = trace_smooth;

            % thresholding
            MAD = mad(trace_smooth, 1);
            thr = 8 * MAD;

            % Peak detection
            [pks, locs] = findpeaks(trace_smooth, 'MinPeakHeight', thr, ...
                'MinPeakProminence', 4*MAD, 'MinPeakDistance', 0.2 * Fs);

            PeakIdx{i,1} = i;
            PeakIdx{i,2} = pks;
            PeakIdx{i,3} = locs;

            signal_mean = mean(pks, 'omitnan');
            SNR = signal_mean / MAD;
            signal_skew = skewness(trace_smooth, 0);

            Metric(i,:) = [i-1, length(pks), SNR, signal_skew];

            % if mod(i,10) == 0
            %     fprintf('Calculating %d / %d\n', i, nCells);
            % end

            %
            % figure
            % subplot(211)
            % hold on
            % plot(full_ts, interp_trace)
            % plot(full_ts, baseline_drift)
            %  plot([0 max(full_ts)], [0 0], 'k--')
            % subplot(212)
            % hold on
            % plot(full_ts, trace_detrended)
            % plot([0 max(full_ts)], [0 0], 'k--')

            % figure
            % figure
            % hold on
            % plot(full_ts, trace_smooth)
            % plot([0 max(full_ts)], [0 0], 'k--')
            % plot([0 max(full_ts)], [thr thr], 'k--')
            % plot(full_ts(locs), pks, 'kv')

            % figure
            % nbins = 100;
            % h = histogram(trace_smooth, nbins);

        end

        Metric(:,1) = Metric(:,1) + 1;

        %% remove low skew ICs
        a = find(Metric(:,4)>1);
        AcceptedMetric_skew = Metric(a,:);


        %% check temporal correlation
        close all

        Pearson_R_all = cell(nCells,1);
        X_R_all = cell(nCells,1);
        idx_close_all = cell(nCells,1);

        for i = 1:nCells
            %%

            % figure
            % plot(full_ts, trace_smooth_all(i,:))

            D_temp = D;
            idx_close = find(D_temp(:,i));

            n = numel(idx_close);
            Pearson_R = nan(n, 1);
            X_R = nan(n, 1);

            A = trace_smooth_all(i,:);
            for kk = 1:n
                k = idx_close(kk);

                B = trace_smooth_all(k,:);

                r = corr(A(:), B(:)); 
                Pearson_R(kk,1) = r;
            end

            Pearson_R_all{i,1} = Pearson_R;
            idx_close_all{i,1} = idx_close;

            % if mod(i,10) == 0
            %     fprintf('Corr Calculating %d / %d\n', i, nCells);
            % end
        end

        %% Remove highly correlated ICs as over-separated cells

        MaxR = cell2mat( cellfun(@(x)max(x), Pearson_R_all, 'UniformOutput', false) );

        R_thr = 0.7;
        Ridx = find(MaxR > R_thr);

        Metric_temp = Metric;
        a = find(Metric_temp(:,4)<=1);
        Metric_temp(a,:) = NaN;

        for i = 1:length(Ridx)
            idx = Ridx(i);

            CloseCellID = idx_close_all{idx};
            R = Pearson_R_all{idx};

            CompareID = [idx, CloseCellID(find(R > R_thr))']; % find ICind to be removed

            SNR = Metric_temp(CompareID, 3);
            [M, I] = max(SNR);
            RemoveID = CompareID; RemoveID(I) =[];

            Metric_temp(RemoveID,:) = NaN;
        end

        AcceptedMetric = Metric_temp;
        rows_with_nan = any(isnan(AcceptedMetric), 2); % Find rows containing any NaN
        AcceptedMetric = AcceptedMetric(~rows_with_nan, :);

        %% show image after removing cells
        close all

        s = strcat('CalB', mice_str_Folder{mice_ind}, " Day", num2str(Day));
        
        fig_originalImage = figure;
        A = ( image_data./prctile(image_data(:), 95) ) * 0.6;
        imshow(A);
        %get(gcf,'Position')
        pos = [845   265   988   662];
        set(gcf,'Position',pos);
        % title('MP image')
        title(strcat(s, " MP image"))

        % original
        fig_originalDetection = figure;
        A = ( image_data./prctile(image_data(:), 95) ) * 0.6;
        I3 = labeloverlay(A, CellImage, 'Colormap','jet','Transparency',0.4);
        imshow(I3); hold on;
        % plot all boundaries
        for kk = 1:nCells
            if ~isempty(Boundaries{kk})
                plot(Boundaries{kk}(:,2), Boundaries{kk}(:,1), 'r', 'LineWidth',0.5);
            end
        end
        set(gcf,'Position',pos);
        title(strcat(s, " ICs original"))


        % remove onle low skewness ICs
        acceptedIDskew = AcceptedMetric_skew(:,1);
        mask = ismember(CellImage, acceptedIDskew);
        CellImage_acceptedCellskew = CellImage .* mask;

        fig_Detection_skew = figure;
        A = ( image_data./prctile(image_data(:), 95) ) * 0.6;
        I3 = labeloverlay(A , CellImage_acceptedCellskew, 'Colormap','jet','Transparency',0.4);
        imshow(I3)
        hold on
        nCells_accepted = length(acceptedIDskew);
        Boundaries_skew = cell(nCells_accepted ,1);
        for kk = 1:nCells_accepted
            img = CellImage_acceptedCellskew == acceptedIDskew(kk);
            CellMask = func_findboundary(img);
            B = bwboundaries(CellMask);
            if ~isempty(B)
                Boundaries_skew{kk} = B{1}; % cell boundaries
                plot(Boundaries_skew{kk}(:,2), Boundaries_skew{kk}(:,1), 'r', 'LineWidth',0.5);
            end
        end
        set(gcf,'Position',pos);
        title(strcat(s, " ICs removed low-skewness traces"))


        % remove high correlated ICs
        acceptedIDs = AcceptedMetric(:,1);
        mask = ismember(CellImage, acceptedIDs);
        CellImage_acceptedCells = CellImage .* mask;

        fig_Detection_skew_corr = figure;
        A = ( image_data./prctile(image_data(:), 95) ) * 0.6;
        I3 = labeloverlay(A , CellImage_acceptedCells, 'Colormap','jet','Transparency',0.4);
        imshow(I3)
        hold on
        nCells_accepted = length(acceptedIDs);
        Boundaries_cleaned = cell(nCells_accepted ,1);
        for kk = 1:nCells_accepted
            img = CellImage_acceptedCells == acceptedIDs(kk);
            CellMask = func_findboundary(img);
            B = bwboundaries(CellMask);
            if ~isempty(B)
                Boundaries_cleaned{kk} = B{1}; % cell boundaries
                plot(Boundaries_cleaned{kk}(:,2), Boundaries_cleaned{kk}(:,1), 'r', 'LineWidth',0.5);
            end
        end
        set(gcf,'Position',pos);
        title(strcat(s, " ICs removed low-skewness and highly-synchronyzed traces"))



        % % temp
        % % original
        % fig_originalDetection = figure;
        % A = ( image_data./prctile(image_data(:), 95) ) * 0.6;
        % imshow(A); hold on;
        % % I3 = labeloverlay(A , CellImage, 'Colormap','jet','Transparency',0.4);
        % % imshow(I3); hold on;
        % % plot all boundaries
        % for kk = [22 132 187 575] %1:nCells
        %     if ~isempty(Boundaries{kk})
        %         plot(Boundaries{kk}(:,2), Boundaries{kk}(:,1), 'r', 'LineWidth',1);
        %     end
        % end
        % set(gcf,'Position',pos);
        % % title('ICs original')

        %% save results

        PeakIdx_accepted = PeakIdx(AcceptedMetric(:,1),:);
        interp_traces_accepted = interp_traces(:, AcceptedMetric(:,1));


        currentFolder = pwd;
        s_savefolder = 'results_Validation';
        [ status, msg ] = mkdir(s_savefolder);
        cd(fullfile(pwd, s_savefolder))

        s = strcat('CellTraces_Validated_', mice_str{mice_ind}, '_d', num2str(Day));
        save(s, 'PeakIdx_accepted', 'AcceptedMetric', 'full_ts', 'interp_traces_accepted', 'image_data',...
            'CellImage_acceptedCells', 'Boundaries_cleaned',...
            '-v7.3');


        ss = strcat(s, '_fig_originalImage');
        print(fig_originalImage, ss, '-dmeta')

        ss = strcat(s, '_fig_originalDetection');
        print(fig_originalDetection, ss, '-dmeta')

        ss = strcat(s, '_fig_Detection_skew');
        print(fig_Detection_skew, ss, '-dmeta')

        ss = strcat(s, '_fig_Detection_skew_corr');
        print(fig_Detection_skew_corr, ss, '-dmeta')


        cd(currentFolder)

    end
end



%% functions

function CellImage = func_findboundary(image)

image(image < 0) = 0;
image = image / sum(image(:));

thrimg = 5 * std(image(:));

image(image < thrimg) = 0;
% image(image >= thrimg) = 1;

a = image(:);
a(a == 0) = [];
P = prctile(a, 80);
thrimg = P;
image(image < thrimg) = 0;
image(image >= thrimg) = 1;

CellImage = image;
end


%% find file name lists
function fnamelist = filenamelisting(path, str_rule)

listing = dir(path);

ftable = struct2table(listing);
A = ftable.name;

clear('fname')
for i = 1:length(A)
    % if  (isempty( regexp(A{i}, str_rule, 'once')) == 0) && strcmp(A{i}([1 2]), 's ') ==0
    if  (isempty( regexp(A{i}, str_rule, 'once')) == 0)
        % if  (isempty( regexp(A{i}, '-(\d)(\d).jpg', 'once')) == 0) ==1
        fname{i,1} = A{i};
        fname{i,2} = ftable.date{i};
        fname{i,3} = ftable.bytes(i);
        fname{i,4} = ftable.datenum(i);
    end
end

x = find(cellfun('isempty',fname(:,1)));
fname(x,:) = [];
fnamelist = fname;

% disp(fname)
end