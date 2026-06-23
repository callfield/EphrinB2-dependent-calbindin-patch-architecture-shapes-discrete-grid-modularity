
function [ Min_meanErrDist Min_maxErrDist]=Fn_run_CorrDecode_ncn_GridMod12(s,t,WorE,Dir,Mod,sort_R2ID,NAME)

             DIR=Dir{s,t}; 
             load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'Grid_Cells');
             All_CELLSET=[Grid_Cells(Mod{s,t}{1}),Grid_Cells(Mod{s,t}{2})];
             
             R2ID = sort_R2ID{WorE,s,t};
%              NAME="GridMod1&2";
        [ Min_meanErrDist Min_maxErrDist]= Fn_CorrD_NCN_GridMod_zRateMap_230328(DIR,All_CELLSET,R2ID,NAME);
        
end