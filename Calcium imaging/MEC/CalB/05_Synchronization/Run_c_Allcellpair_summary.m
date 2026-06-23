clear all;close all;


addpath(pwd)

addpath("function")
CDir=pwd;
load(fullfile(CDir, "..", "Data.mat"), "iEB2"); 
load(fullfile(CDir, "Data.mat"), "SigfCorr_ratio", "EB2_SigfCorr_ratio");
load(fullfile(CDir, "Island.mat"), "wt_Island", "eb_Island") % manually define island cell; 
% wt_Island{s,t} is 0, 1,2,~ 0: no-island cell, 1: island 1, 2: island 2...
% eb_Island{s,t} is 0, 1  0: no-island cell, 1: island cell

% Posi_SigfCorr_ratio=cell(5,3);
% Nega_SigfCorr_ratio=cell(5,3);
% Posi_SigfCorrMean=cell(5,3);
% Nega_SigfCorrMean=cell(5,3);
% Posi/Nega_SigfCorr_ratio: ratio of sig. corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
% 
% Posi/Nega_SigfCorrMean: mean of sig. positive corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
%
% NETWORK: cell pair matrix for each network
% 1:Grid, 2:nonGrid, 3:IntraModule(if exist), 4:InterModule

mkdir Sum

%% Ratio
% 1: ALL-ALL
% 2: Intra Island(I1-I1, I2-I2, I3-I3)
% 3: non Island
% 4: Trans Island(I1-I2, I1-I3, I2-I3)

ALL=[];Intra=[];nonIsland=[];Trans=[];
i=1;
for s=1:5
    for t=1:5
         disp(t)
        if isempty(wt_Island{s,t})==0

            DATA=SigfCorr_ratio{s,t}(:,:,:,6);%r>0.2
            ISLAND=wt_Island{s,t};
            [ALL{i}, Intra{i}, nonIsland{i}, Trans{i}]=Fn_SubCorr_Island_cellpair(DATA,ISLAND);
            i=i+1;

        end
    end
end


ALLmat = vertcat(ALL{:});  
Intramat = vertcat(Intra{:});  
Transmat = vertcat(Trans{:}); 

writematrix(ALLmat, 'Corr_cellpair.xlsx','Sheet', 'WT_0.2_all',"AutoFitWidth",false)
writematrix(Intramat, 'Corr_cellpair.xlsx','Sheet', 'WT_0.2_intra',"AutoFitWidth",false)
writematrix(Transmat, 'Corr_cellpair.xlsx','Sheet', 'WT_0.2_trans',"AutoFitWidth",false)





ALLeb=[];Intra_eb=[];Top10_eb=[];Trans_eb=[];
i=1;
for s=1:5
    for t=1:3
        if isempty(eb_Island{s,t})==0
            mIslandXY=iEB2{s,t}(find(eb_Island{s,t}==1),1:2);
            tmp=pdist2(mIslandXY,mIslandXY);
            TMPNET=nan(size(tmp,1),size(tmp,2));
            TMPNET(tmp<250/2)=1; % select cell set within 250/2um

            DATA=EB2_SigfCorr_ratio{s,t}(:,:,:,6);%r>0.2
             ISLAND=eb_Island{s,t};
             DATA2=iEB2{s,t};
           
            [ALLeb{i}, Intra_eb{i}, Top10_eb{i}, Trans_eb{i}]= Fn_SubCorr_psuedoIsland_cellpair( DATA, ISLAND, TMPNET, DATA2);
            
            i=i+1;
        end
    end
end




ALLeb_mat = vertcat(ALLeb{:});  
Intra_eb_mat = vertcat(Intra_eb{:});  
Top10_eb_mat = vertcat(Top10_eb{:});  
Trans_eb_mat = vertcat(Trans_eb{:}); 

writematrix(ALLeb_mat, 'Corr_cellpair_allEB2.csv');
writematrix(ALLeb_mat, 'Corr_cellpair.xlsx','Sheet', 'EB2_0.2_all',"AutoFitWidth",false)
writematrix(Intra_eb_mat, 'Corr_cellpair.xlsx','Sheet', 'EB2_0.2_pintra',"AutoFitWidth",false)
writematrix(Trans_eb_mat, 'Corr_cellpair.xlsx','Sheet', 'EB2_0.2_ptrans',"AutoFitWidth",false)
writematrix(Top10_eb_mat, 'Corr_cellpair.xlsx','Sheet', 'EB2_0.2_Top10pintra',"AutoFitWidth",false)






% Example: ALLmat and ALLeb_mat are combined into one vector with group labels.

% Stack data as a vertical vector.

data = [ALLmat; ALLeb_mat];

% Create group labels.
group_num = [ ...
    ones(length(ALLmat), 1); ...
    2 * ones(length(ALLeb_mat), 1)];

% Convert numeric labels to categorical labels.
group = categorical(group_num);
% Violinplot
figure('Position', [100, 100, 400, 500]);
violinplot( group, data);


ylim([0 100])

box off;

% Save as PDF.
output_pdf = 'ViolinPlot_ALLmat_vs_ALLeb_mat.pdf';
exportgraphics(gcf, output_pdf, 'ContentType', 'vector', 'Resolution', 300);

