function [ sort_r2ID]=Fn_R2_Contribution_zRateMap_220328(DIR)

%   CELLSET=Grid_Cells;
%    DIR=wt_Dir{s,t}; 

load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')


XXXXX=split(DIR,"\");
samplename=strcat(char(strrep(XXXXX(end-1),'_','-')), " ",char(strrep(strrep(XXXXX(end),'_',' '),'OF','')))
XXXXX=[];
CDir=pwd;


numCells=size(GSrate_map,1);


%%
tmpGSrate_map=nan(numCells,50,50);
for c=1:numCells
    tmpGSrate_map(c,:,:)=GSrate_map{c};
end

% z scored
zTotalRatemap=nan(numCells,50,50);
    MEAN=mean(tmpGSrate_map,[2,3],"omitnan");
    STD=std(tmpGSrate_map,[],[2 3],"omitnan");
    zTotalRatemap=(tmpGSrate_map-MEAN)./STD;

%  disp('corrmap start')

for BaseY=1:50
    for y=1:50

        tmp=corr(zTotalRatemap(:,:,BaseY), ...
            zTotalRatemap(:,:,y), 'rows','complete');
        for x=1:50
            CORRMAP(:,BaseY,x,y)=tmp(:,x);
        end
    
    end
end
    disp('corrmap end')



sort_r2=cell(50,50);
sort_r2ID=cell(50,50);
for BaseX=1:50

    if mod(BaseX,25) ==1
        BaseX
    end
% tic
    parfor (BaseY=1:50,2)

         [ sort_r2ID{BaseX,BaseY}]=Fn_R2_220323(zTotalRatemap, CORRMAP,BaseX,BaseY)
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
    
    % toc
end

end
