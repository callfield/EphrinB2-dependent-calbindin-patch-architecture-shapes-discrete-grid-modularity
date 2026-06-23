
function [ Min_meanErrDist Min_maxErrDist]=Fn_run_CorrDecode_ncn_allGrid(s,t,WorE,Dir,Mod,sort_r2ID,NAME)

             DIR=Dir{s,t}; 
             load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'Grid_Cells');
             All_CELLSET=Grid_Cells;
             
             R2ID = sort_r2ID{WorE,s,t};
%              NAME="GridMod1&2";
        [ Min_meanErrDist Min_maxErrDist]= Fn_CorrDecode_NormCellNum_GridMod_230324(DIR,All_CELLSET,R2ID,NAME);
        
end