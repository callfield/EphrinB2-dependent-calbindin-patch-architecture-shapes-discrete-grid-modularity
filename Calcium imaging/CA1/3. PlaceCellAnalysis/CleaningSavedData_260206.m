clear
close all

%% load data
% mice_str = {'p13', 'p14','p15', 'p16'};
% mice_str_Folder = {'P13', 'P14','P15', 'P16'};
mice_str = {'n12', 'n14'};
mice_str_Folder = {'N12', 'N14'};

Nanimal = 2;
Ndays = 6;

thr = 6;

L = cell(Nanimal, Ndays);
for Day = 1:6
    for mice_ind = 1:2
        s_loadfolder_base = 'H:\experiments H drive\260121 Ca imaging\Behavior\PlaceFieldRawAnalysis_260204\6MAD';
        path_load_animal = fullfile(s_loadfolder_base, mice_str_Folder{mice_ind});
        path_load_Day = fullfile(path_load_animal, strcat('Day', num2str(Day)));

        baseFolder = 'H:\experiments H drive\260121 Ca imaging\Behavior\PlaceFieldRawAnalysis_260204';
        saveFolder = fullfile(baseFolder, ...
            sprintf('%dMAD', thr), mice_str_Folder{mice_ind}, ...
            sprintf('Day%d', Day));
        saveFolder = fullfile(saveFolder, 'Basics_260204');

        s_loadmat = strcat('SpaceAnalysis_', mice_str{mice_ind}, '_d', num2str(Day), '.mat');
        ss = fullfile(saveFolder, s_loadmat);


        if exist(ss, 'file') == 0
            fprintf(strcat("File not found: ", 'Mouse', mice_str_Folder{mice_ind}, " Day%d\n"), Day);
            continue;
        end

        fprintf(sprintf('loading %s\n', s_loadmat));
        L = load(ss);

        fprintf('start cleaning\n');


        %% deleat unusing fields in structure

        LData_ori = L;
        LData = L;
        nCell = length(LData.SpInfo_cell);

        T = structSizeTable(LData);


        LData.Detections.AnimalSpeed = LData.RateMaps.SpkDetection.SpeedResults.smoothedSpeed;
        LData.Detections.Dur_Running = LData.RateMaps.Adapt.Dur_Running;
        LData.Detections.numEvent = LData.RateMaps.SpkDetection.numEvent;

        % RateMap images will not be used in later analysis
        LData = rmfield(LData, 'RateMaps');

        % New FieldAnalysis was saved in other files
        LData = rmfield(LData, 'FieldAnalysis_cell');

        % Remove images since their sizes are large
        FieldName = {'XCor_sp', 'rateMapUniform_Gaussian_sp_halve_all', 'rateMapUniform_sp_halve_all', 'occMapUniform_Gaussian_sp_halve_all'};
        for i_Cell = 1:nCell
            for i = 1:length(FieldName)
                LData.HalveCompareField_stability_cell{i_Cell}      = rmfield(LData.HalveCompareField_stability_cell{i_Cell}, FieldName{i});
                LData.HalveCompareField_reliability_cell{i_Cell}    = rmfield(LData.HalveCompareField_reliability_cell{i_Cell}, FieldName{i});
            end
        end

        %% Remove images in shuffle stability
        LData_bak = LData;

        FieldName = {'XCor_sp','rateMapUniform_sp_halve_all', 'rateMapUniform_Gaussian_sp_halve_all', 'occMapUniform_Gaussian_sp_halve_all'};
        for i_Cell = 1:nCell
            for i = 1:length(FieldName)
                LData.Halves_shuf_stab_cell{i_Cell}.Field.real     = rmfield(LData.Halves_shuf_stab_cell{i_Cell}.Field.real, FieldName{i});
                LData.Halves_shuf_reliab_cell{i_Cell}.Field.real   = rmfield(LData.Halves_shuf_reliab_cell{i_Cell}.Field.real, FieldName{i});
            end
        end
        i_Cell = 1;
        TT2 = structSizeTable(LData.Halves_shuf_stab_cell{i_Cell});
        T2 = structSizeTable(LData);

        % S = whos(A);
        % field = 'RateMaps';
        % A = rmfield(A,field);

        BasicAnalysis = LData;


        %% Clean Field ReAnalized data

        s_loadmat_field = strcat('FieldAnalysisRedo_', mice_str{mice_ind}, '_d', num2str(Day), '.mat');
        ss = fullfile(path_load_Day, 'FieldRedo_260204_Ratio', s_loadmat_field);

        fprintf(sprintf('loading %s\n', s_loadmat_field));
        F_ori = load(ss);

        % fprintf('start cleaning\n');

        %%

        F = F_ori;
        for i_Cell = 1:nCell
            F_temp = F.FRD_FieldAnalysis_Redo_cell{i_Cell};
            FieldName = {'original_RateMap', 'thres_RateMap', 'watershed_RateMap', 'AllThrImages' };
            F_temp = removevars(F_temp, FieldName);
            F.FRD_FieldAnalysis_Redo_cell{i_Cell} = F_temp;
        end

        for i_Cell = 1:nCell
            F_temp = F.FRD_ACGanalysis_cell{i_Cell};
            % T = varSizeTable(F_temp);
            % FieldName = {'Autoc'};
            % F_temp = removevars(F_temp, FieldName);
            F_temp   = rmfield(F_temp, 'Autoc');
            F.FRD_ACGanalysis_cell{i_Cell} = F_temp;
        end

        % i_Cell = 5;
        % F_temp = F.ACGanalysis_cell{i_Cell};
        % T = varSizeTable(F_temp);

        FieldAnalysisRedo = F;

        %% bootstrap

        s_loadmat_boot = strcat('RateMapBootstrap_Shuffle', mice_str{mice_ind}, '_d', num2str(Day), '.mat');
        ss = fullfile(path_load_Day, 'RateMapBootstrap_shuffle_260204', s_loadmat_boot);

        fprintf(sprintf('loading %s\n', s_loadmat_boot));
        BOOT = load(ss);

        % fprintf('start cleaning\n');
        % T = structSizeTable(BOOT);

        %% topology

        Dir = 'H:\experiments H drive\251203 Ca imaging\code_251211\CellValidation\Control\results_Validation';
        s_loadmat_topology = strcat('CellTraces_Validated_', mice_str{mice_ind}, '_d', num2str(Day), '.mat');
        ss = fullfile(Dir, s_loadmat_topology);

        fprintf(sprintf('loading %s\n', s_loadmat_topology));
        Topo = load(ss);

        % T = structSizeTable(Topo);
        Topo = rmfield(Topo, 'interp_traces_accepted');


        %% save centroides

        cleanedBoundaries = cellfun(@(x) unique(x(~any(isnan(x),2),:),'rows','stable'), ...
            Topo.Boundaries_cleaned, 'UniformOutput', false);
        polyin = cellfun(@(x) polyshape(x, 'Simplify', true), cleanedBoundaries);
        [y, x] = centroid(polyin);
        CellCenter = [x, y];
        Topo.CellCenter = CellCenter;
        Topo.umperpx = 3.2;


        %% save cleaned data

        Folder = 'H:\experiments H drive\251203 Ca imaging\code_251211\Place Field analysis\statistics\CleanedAnalysisData';
        OutDir = fullfile(Folder, 'CleanedAnalysisData_AdpBin_ThrPeakRatio');
        [status,msg,msgID] = mkdir(OutDir);

        s_savemat = strcat('CleanedAnalysisData_', mice_str{mice_ind}, '_d', num2str(Day));
        ss = fullfile(OutDir, s_savemat);

        save(ss, 'BasicAnalysis', 'FieldAnalysisRedo', 'BOOT', 'Topo', ...
            '-v7.3');


    end
end





%% function
%% investigate size of saved field in struct
function T = structSizeTable(S, prefix)
if nargin < 2
    prefix = '';
end
T = table();

fields = fieldnames(S);
for i = 1:numel(fields)
    f = fields{i};
    value = S.(f);
    fullName = [prefix '.' f];
    if startsWith(fullName, '.')
        fullName = fullName(2:end);
    end

    if isstruct(value)
        subT = structSizeTable(value, fullName);
        T = [T; subT]; %#ok<AGROW>
    else
        % get size
        tmp = value; %#ok<NASGU>
        info = whos('tmp');
        newRow = {fullName, class(value), size(value), info.bytes, info.bytes/1024/1024};
        T = [T; cell2table(newRow, 'VariableNames', ...
            {'Field','Class','Size','Bytes','MB'})]; %#ok<AGROW>
    end
end
end



%%
function T = varSizeTable(X, prefix)
if nargin < 2
    prefix = inputname(1);
    if isempty(prefix), prefix = 'ans'; end
end
rows = {};

    function pushRow(fieldName, className, sz, bytes)
        rows(end+1,:) = {fieldName, className, mat2str(sz), bytes, bytes/1024/1024};
    end

% ---- struct ----
if isstruct(X)
    if isempty(X)
        tmp = X; info = whos('tmp'); %#ok<NASGU>
        pushRow(prefix, class(X), size(X), info.bytes);
    elseif numel(X) > 1
        for idx = 1:numel(X)
            subName = sprintf('%s(%d)', prefix, idx);
            subT = varSizeTable(X(idx), subName);
            rows = [rows; table2cell(subT)]; %#ok<AGROW>
        end
    else
        flds = fieldnames(X);
        for k = 1:numel(flds)
            f = flds{k};
            if isempty(X(1).(f))
                val = X(1).(f);
                tmp = val; info = whos('tmp'); %#ok<NASGU>
                pushRow(sprintf('%s.%s',prefix,f), class(val), size(val), info.bytes);
            else
                val = X(1).(f);
                subName = sprintf('%s.%s', prefix, f);
                subT = varSizeTable(val, subName);
                rows = [rows; table2cell(subT)]; %#ok<AGROW>
            end
        end
    end

    % ---- table ----
elseif istable(X)
    if isempty(X)
        tmp = X; info = whos('tmp'); %#ok<NASGU>
        pushRow(prefix, class(X), size(X), info.bytes);
    else
        vars = X.Properties.VariableNames;
        for k = 1:numel(vars)
            vname = vars{k};
            val = X.(vname);
            subName = sprintf('%s.%s', prefix, vname);
            subT = varSizeTable(val, subName);
            rows = [rows; table2cell(subT)]; %#ok<AGROW>
        end
    end

    % ---- cell ----
elseif iscell(X)
    if isempty(X)
        tmp = X; info = whos('tmp'); %#ok<NASGU>
        pushRow(prefix, class(X), size(X), info.bytes);
    else
        for k = 1:numel(X)
            val = X{k};
            subName = sprintf('%s{%d}', prefix, k);
            subT = varSizeTable(val, subName);
            rows = [rows; table2cell(subT)]; %#ok<AGROW>
        end
    end

    % ---- others ----
else
    tmp = X; info = whos('tmp'); %#ok<NASGU>
    pushRow(prefix, class(X), size(X), info.bytes);
end

if isempty(rows)
    T = table([],[],[],[],[],'VariableNames',{'Field','Class','Size','Bytes','MB'});
else
    T = cell2table(rows, 'VariableNames', {'Field','Class','Size','Bytes','MB'});
end
end


