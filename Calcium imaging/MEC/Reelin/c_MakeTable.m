close all; clear all;




% MakeDirFile.m で前もってgroup1_Dir, group2_Dirを作成しなおす
load('Data.mat')
load('04_Head_Direction\Data_HD.mat')
load('02_Spatial_cell\Spatial_Info.mat')
load('08_GridModule\GridMod.mat');

GROUP1_T = cell(5,3);   

for s = 1:5
    for t = 1:3

        DIR = group1_Dir{s,t};

        %==========================================================
        % Load 基本データ
        %==========================================================
        load(fullfile(DIR,'ST_dF_grid_aut_data.mat'), ...
            'Original_Cell_ID','Grid_Cells','GSrate_map', ...
            'B_Score_v2','Border_cells_v2','RecStart','sTrk')

        nCell = numel(Original_Cell_ID);

        %==========================================================
        % Table 初期化（ここは一時変数 TT を使う）
        %==========================================================
        TT = table;
        TT.AnimalName          = repmat(SampleName(1,s), nCell, 1);
        TT.Trial               = repmat(strcat("Trial",num2str(t)), nCell, 1);
        TT.Cell_ID             = (1:nCell).';
        TT.Original_Cell_ID    = Original_Cell_ID(:);
        TT.DV_position_um       = GROUP1{s,t}(:,1); 
        TT.ML_position_um       = GROUP1{s,t}(:,2); 

        %==========================================================
        % Grid / HD / Spatial / Border スコア
        %==========================================================
        TT.Grid_Score       = GROUP1{s,t}(:,13);
        TT.Grid_scale       = GROUP1{s,t}(:,3);
        TT.Grid_orientaion  = GROUP1{s,t}(:,4);
        TT.Grid_field       = GROUP1{s,t}(:,5);
        TT.Grid_cell        = GROUP1{s,t}(:,3) > 0;

        TT.Border_Score     = B_Score_v2;
        TT.HD_Score         = group1_HD_Score{{s,t};
% 1 : HD score（Conditioned entropy from Fireing rateの(6degree Bin )）
% 2 : Watson's U²-test p-value, （ Watson's U²-test : a nonparametric permutation test based on Watson's U2 statistic for 2-sample testing of circular data (angles).
% 3 : U2 Score
% 4: RL test  p-value, (Rayleigh test（null: uniformity, alt: unimodal directionality）)
% 5: RL test z statistics
% 6: Angular Standard Deviation
        TT.Spatial_Score    = Spatial_Info_GROUP1{s,t}(:,3);

        %==========================================================
        % Grid module（0/1/2/3）
        %==========================================================
        TT.Grid_module = zeros(nCell,1);

        if  isempty(gp1_Mod{s,t}) == 0;
            Mod1 = Grid_Cells(gp1_Mod{s,t}{1});
            Mod2 = Grid_Cells(gp1_Mod{s,t}{2});
            Othr = setdiff(Grid_Cells, [Mod1; Mod2]);
    
            TT.Grid_module(Mod1) = 1;
            TT.Grid_module(Mod2) = 2;
            TT.Grid_module(Othr) = 3;
        else
            TT.Grid_module(Grid_Cells) = 3;
        end

        %==========================================================
        % Rate Map（リンク）
        %==========================================================
        rm_path = fullfile(DIR,'GSrate_map.mat');
        if ~isfile(rm_path)
            save(rm_path, 'GSrate_map','-v7.3');
        end
        mf_rm = matfile(rm_path,'Writable',false);

        TT.RM_idx    = (1:nCell).';
        TT.RM_source = repmat({mf_rm}, nCell, 1);

        %==========================================================
        % dF（リンク）
        %==========================================================
        dF = readmatrix(fullfile(DIR,'ST_PCI_noDup_dF.csv'));

        df_path = fullfile(DIR,'dF_data.mat');
        if ~isfile(df_path)
            save(df_path,'dF','-v7.3');
        end
        mf_df = matfile(df_path,'Writable',false);

        TT.dF_col    = (1:nCell).';
        TT.dF_source = repmat({mf_df}, nCell, 1);

        %==========================================================
        % Track & Timestamp（リンク）
        %==========================================================
        Trk_withTimeStamp = makeTimeStamp(DIR);

        ts_path = fullfile(DIR,'Trk_wTS.mat');
        if ~isfile(ts_path)
            save(ts_path,'Trk_withTimeStamp','-v7.3');
        end

        mf_ts = matfile(ts_path,'Writable',false);

        TT.Trk_idx    = (1:nCell).';
        TT.Trk_source = repmat({mf_ts}, nCell, 1);


        GROUP1_T{s,t} = TT;



    end
end

save('Group1_CellTables.mat', 'GROUP1_T', '-v7.3');

GROUP2_T = cell(5,3);  

for s = 1:7
    for t = 1:3

    
            DIR = group2_Dir{s,t};
    
            %==========================================================
            % Load 基本データ
            %==========================================================
            load(fullfile(DIR,'ST_dF_grid_aut_data.mat'), ...
                'Original_Cell_ID','Grid_Cells','GSrate_map', ...
                'B_Score_v2','Border_cells_v2','RecStart','sTrk')
    
            nCell = numel(Original_Cell_ID);
    
    
            %==========================================================
            % Table 初期化
            %==========================================================
            T = table;
            T.AnimalName          = repmat(SampleName(2,s), nCell, 1);
            T.Trial               = repmat(strcat("Trial",num2str(t)), nCell, 1);
            T.Cell_ID             = (1:nCell).';
            T.Original_Cell_ID    = Original_Cell_ID(:);
            T.DV_position_um       = GROUP2{s,t}(:,1); 
            T.ML_position_um       = GROUP2{s,t}(:,2); 
    
            %==========================================================
            % Grid / HD / Spatial / Border スコア
            %==========================================================
            T.Grid_Score       = GROUP2{s,t}(:,13);
            T.Grid_scale       = GROUP2{s,t}(:,3);
            T.Grid_orientaion  = GROUP2{s,t}(:,4);
            T.Grid_field       = GROUP2{s,t}(:,5);
            T.Grid_cell        = GROUP2{s,t}(:,3) > 0;
    
            T.Border_Score     = B_Score_v2;
            T.HD_Score         = group2_HD_Score{s,t};
% 1 : HD score（Conditioned entropy from Fireing rateの(6degree Bin )）
% 2 : Watson's U²-test p-value, （ Watson's U²-test : a nonparametric permutation test based on Watson's U2 statistic for 2-sample testing of circular data (angles).
% 3 : U2 Score
% 4: RL test  p-value, (Rayleigh test（null: uniformity, alt: unimodal directionality）)
% 5: RL test z statistics
% 6: Angular Standard Deviation
            T.Spatial_Score    = Spatial_Info_GROUP2{s,t}(:,3);
    
    
            %==========================================================
            % Grid module（0/1/2/3）
            %==========================================================
            T.Grid_module = zeros(nCell,1);
            if  isempty(gp2_HD_Score{s,t}) == 0;
                Mod1 = Grid_Cells(gp2_HD_Score{s,t}{1});
                Mod2 = Grid_Cells(gp2_HD_Score{s,t}{2});
                Othr = setdiff(Grid_Cells, [Mod1; Mod2]);
        
                T.Grid_module(Mod1) = 1;
                T.Grid_module(Mod2) = 2;
                T.Grid_module(Othr) = 3;
            else
                T.Grid_module(Grid_Cells) = 3;
            end
    
            %==========================================================
            % Rate Map（リンク）
            %==========================================================
            rm_path = fullfile(DIR,'GSrate_map.mat');
            if ~isfile(rm_path)
                save(rm_path, 'GSrate_map','-v7.3');
            end
            mf_rm = matfile(rm_path,'Writable',false);
    
            T.RM_idx    = (1:nCell).';
            T.RM_source = repmat({mf_rm}, nCell, 1);
    
    
            %==========================================================
            % dF（リンク）
            %==========================================================
            dF = readmatrix(fullfile(DIR,'ST_PCI_noDup_dF.csv'));
    
            df_path = fullfile(DIR,'dF_data.mat');
            if ~isfile(df_path)
                save(df_path,'dF','-v7.3');
            end
            mf_df = matfile(df_path,'Writable',false);
    
            T.dF_col    = (1:nCell).';
            T.dF_source = repmat({mf_df}, nCell, 1);
    
    
            %==========================================================
            % Track & Timestamp（リンク）
            %==========================================================
            Trk_withTimeStamp = makeTimeStamp(DIR);
    
            ts_path = fullfile(DIR,'Trk_wTS.mat');
            if ~isfile(ts_path)
                save(ts_path,'Trk_withTimeStamp','-v7.3');
            end
    
            mf_ts = matfile(ts_path,'Writable',false);
    
            T.Trk_idx    = (1:nCell).';
            T.Trk_source = repmat({mf_ts}, nCell, 1);
    
            GROUP2_T{s,t} = T;


    end
end

save('Group2_CellTables.mat', 'GROUP2_T', '-v7.3');







