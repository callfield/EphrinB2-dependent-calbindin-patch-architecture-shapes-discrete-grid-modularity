function [ sort_r2ID]=Fn_CorrDecode_Contribution_220323(DIR)
% すべて細胞のratematを重ね合わせたものを作成（Maxの正確性があるはず）
% 重ね合わせたrateマップを用いて各spatial　binに対して、corrmapを作成。
% すべての細胞を重ね合わせたcorrmapにたいして、それぞれの細胞のR2貢献度を作成する
%　R2貢献度の高い順に細胞をソートする
%   CELLSET=Grid_Cells;
%    DIR=wt_Dir{s,t}; 

load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')


XXXXX=split(DIR,"\");
samplename=strcat(char(strrep(XXXXX(end-1),'_','-')), " ",char(strrep(strrep(XXXXX(end),'_',' '),'OF','')))
XXXXX=[];
CDir=pwd;


numCells=size(GSrate_map,1);


%%
% すべて細胞のratemat（Z-score化する）を重ね合わせたものを作成（Maxの正確性があるはず）
TotalRatemap=nan(numCells,50,50);
for c=1:numCells
    TotalRatemap(c,:,:)=GSrate_map{c};
end



% 重ね合わせたrateマップを用いて各spatial　binに対して、corrmapを作成。
for BaseY=1:50
    for y=1:50

        tmp=corr(TotalRatemap(:,:,BaseY), ...
            TotalRatemap(:,:,y), 'rows','complete');
        for x=1:50
            CORRMAP(:,BaseY,x,y)=tmp(:,x);% corrleartion of rate map on 
        end

    end
end
    disp('corrmap end')
% tmp=[];
% tmp(1:50,1:50)=CORRMAP(10,40,:,:);
% imshow(tmp)

sort_r2=cell(50,50);
sort_r2ID=cell(50,50);
for BaseX=1:50

    if mod(BaseX,25) ==1
        BaseX
    end
tic
    parfor (BaseY=1:50,2)
% すべての細胞を重ね合わせたcorrmapにたいして、各binにおいてそれぞれの細胞のR2貢献度を作成する
         [ sort_r2ID{BaseX,BaseY}]=Fn_R2_220323(TotalRatemap, CORRMAP,BaseX,BaseY)
%         ttt=reshape(CORRMAP(BaseX,BaseY,:,:),50,50);
%         Y=reshape(ttt,[],50*50);
% 
%         X=reshape(TotalRatemap(:,:,:),[],50*50);
%         tmp_r2=[];
%         for c=1:numCells
%             mdl = fitlm(X(c,:),Y);
%             tmp_r2(c)=mdl.Rsquared.Ordinary;
%         end
%         [ sort_r2ID{BaseX,BaseY}]=sort(tmp_r2,'descend');

    end
    
    toc
end

end
