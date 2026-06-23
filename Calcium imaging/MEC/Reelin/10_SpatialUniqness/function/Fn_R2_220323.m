function [ sort_r2ID]=Fn_R2_220323(TotalRatemap, CORRMAP,BaseX,BaseY)
        numCells=size(TotalRatemap,1);
        ttt=reshape(CORRMAP(BaseX,BaseY,:,:),50,50);
        Y=reshape(ttt,[],50*50);

        X=reshape(TotalRatemap(:,:,:),[],50*50);
        tmp_r2=[];
        for c=1:numCells
            mdl = fitlm(X(c,:),Y);
            tmp_r2(c)=mdl.Rsquared.Ordinary;
        end
        [sort_r2, sort_r2ID]=sort(tmp_r2,'descend');

end