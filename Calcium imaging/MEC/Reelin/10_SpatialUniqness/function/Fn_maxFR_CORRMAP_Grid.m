function [CORRMAP CORRMAP_Pval]=Fn_maxFR_CORRMAP(num_SHUFFLE,BaseX,BaseY,num_SelectCell,AnaCellset,tmpGSrate_map)
% 各spatial　binに対して,各spatial　bin と他のbinのFireing rate 差分をみる
% tmpGSrate_map=Z_tmpGSrate_map;
% diffTotalRatemap=abs(tmpGSrate_map(:,:,:)-tmpGSrate_map(:,BaseX,BaseY));

%Top active cell on [BaseX BaseY]
[dFR, sortCELLSET]=sort(tmpGSrate_map(AnaCellset,BaseX,BaseY),"descend");
BaseTopCELLSET=sortCELLSET(1:num_SelectCell);


            % tmp_cellset=SHUFFLE + num_SHUFFLE*(SHUFFLE2-1);
            % CELLSET=NC_CELLSET{BaseX,BaseY}(SHUFFLE,:);
CORRMAP=nan(num_SHUFFLE,50,50);
CORRMAP_Pval=nan(num_SHUFFLE,50,50);
        for x=1:50
            for y=1:50
                for SHUFFLE=1:num_SHUFFLE
                    tmpCOR=nan(2,1);
                    tmpPval=nan(2,1);

                   %  Top active cell on [x y] 
                    [dFR, sortCELLSET]=sort(tmpGSrate_map(AnaCellset,x,y),"descend");
                    xyTopCELLSET=sortCELLSET(1:num_SelectCell);

% add Top10 cell on [BaseX BaseY]
                    [tmpCOR(1) tmpPval(1)]=corr(tmpGSrate_map(BaseTopCELLSET,BaseX,BaseY),tmpGSrate_map(BaseTopCELLSET,x,y), 'rows','pairwise');
% add Top10 cell on [x y]
                    [tmpCOR(2) tmpPval(2)]=corr(tmpGSrate_map(xyTopCELLSET,BaseX,BaseY),tmpGSrate_map(xyTopCELLSET,x,y), 'rows','pairwise');


                   CORRMAP(SHUFFLE,x,y)=min(tmpCOR(:),[],"all","omitmissing");
                   if isnan(CORRMAP(SHUFFLE,x,y))==0
                        CORRMAP_Pval(SHUFFLE,x,y)=unique(tmpPval(find(tmpCOR==CORRMAP(SHUFFLE,x,y))));
                   end


                   
                end
            end
        end
end

