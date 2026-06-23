function saveAngle(DIR, pathCSV, preFix)
    % Move to the specified directory
    cd(DIR);

    % Load recording metadata
    load("ST_dF_grid_aut_data.mat", "RecStart");
    vNum = numel(RecStart);


    % Get the current directory information
    [pathParent, currentDirName] = fileparts(pwd);
    [~, parentDirName] = fileparts(pathParent);


    % Containers for angle data
    HeadAngle = [];
    BodyAngle = [];

    % Process each video (trial)
    newV1 = 1; newV2 = 1;
    for v = 1:vNum
        if isempty(RecStart{v})
            % Skip if this trial has no valid start info
        else
            % Compose the CSV file path
            DLC_file = fullfile(pathCSV, sprintf("%s_%s_Trial_%d%s", ...
                               parentDirName, currentDirName, v, preFix));

            % Read CSV data
            DLC_tmp = readmatrix(DLC_file);

            % Compute angles and append to buffers
            % (Head: columns 2-3 and 5-6; Body: columns 8-9 and 11-12)
            [HeadAngle, newV1] = extractAngle(DLC_tmp, RecStart{v}, newV1, 2, 3, 5, 6, HeadAngle);
            [BodyAngle, newV2] = extractAngle(DLC_tmp, RecStart{v}, newV2, 8, 9, 11, 12, BodyAngle);
        end
    end

    % Save angles to file
    save('Angle.mat', 'HeadAngle', 'BodyAngle');
end
