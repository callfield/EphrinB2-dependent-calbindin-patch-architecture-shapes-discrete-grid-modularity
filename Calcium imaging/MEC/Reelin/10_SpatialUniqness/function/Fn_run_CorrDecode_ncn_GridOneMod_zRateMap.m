function [ Min_meanErrDist Min_maxErrDist]=Fn_run_CorrDecode_ncn_GridOneMod(s,t,WorE,Dir,Mod,sort_R2ID,NAME)

             DIR=Dir{s,t}; 
             load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'Grid_Cells');
             All_CELLSET=Grid_Cells(Mod{s,t}{1});         
             R2ID = sort_R2ID{WorE,s,t};
%              NAME="GridMod1&2";
        [ Min_meanErrDist1 Min_maxErrDist1]= Fn_CorrD_NCN_GridMod_zRateMap_230328(DIR,All_CELLSET,R2ID,NAME);
       
        
              All_CELLSET=Grid_Cells(Mod{s,t}{2});         
        [ Min_meanErrDist2 Min_maxErrDist2]= Fn_CorrD_NCN_GridMod_zRateMap_230328(DIR,All_CELLSET,R2ID,NAME);
 
        
        tmpMin_meanErrDist(1,1:50,1:50)=Min_meanErrDist1;
        tmpMin_meanErrDist(2,1:50,1:50)=Min_meanErrDist2;
        
        tmpMin_maxErrDist(1,1:50,1:50)=Min_maxErrDist1;
        tmpMin_maxErrDist(2,1:50,1:50)=Min_maxErrDist2;
        
        
        % minumum error
Min_meanErrDist=reshape(min(tmpMin_meanErrDist,[],1,"omitnan"),[],50);
Min_maxErrDist=reshape(min(tmpMin_maxErrDist,[],1,"omitnan"),[],50);

%%

x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))
x=[];
num_AnaCell=10; % number of one cell set
CDir=pwd;


imagesc(Min_meanErrDist,'AlphaData',~isnan(Min_meanErrDist));% nan to transparent
set(gca,'color', [1 1 1]); % set backgound color
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet); 
cb=colorbar;

ylabel(cb,"Min average error distance (cm)",'FontSize',11);
title({samplename;...
strcat("Min average error (",num2str(num_AnaCell), ...
     " cell from ",NAME,")")},...
    'FontSize',14);

exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\GridModule\Min_meanError\",NAME," Min_meanError " ,samplename,".jpg"))
        clf
    
        
imagesc(Min_maxErrDist,'AlphaData',~isnan(Min_maxErrDist));% nan to transparent
set(gca,'color', [1 1 1]); % set backgound color
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet); 
cb=colorbar;

ylabel(cb,"Min largest error distance (cm)",'FontSize',11);
title({samplename;...
strcat("Min largest error(",num2str(num_AnaCell), ...
     " cell from ",NAME,")")},...
    'FontSize',14);

exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\GridModule\Min_maxError\",NAME," Min_maxError " ,samplename,".jpg"))
        clf
end