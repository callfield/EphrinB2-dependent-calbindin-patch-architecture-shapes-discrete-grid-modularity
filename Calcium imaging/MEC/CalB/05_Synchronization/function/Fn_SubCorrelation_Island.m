function [SigfCorr_ratio, Posi_SigfCorrMean, Nega_SigfCorrMean]=...
    Fn_SubCorrelation_Island(CDir,DIR,samplename,s,t,Island,sWIN)
caFr=10;
save(strcat(DIR,"\ST_dF_grid_aut_data.mat"),"caFr","-append");
%%
% SigfCorr_ratio: ratio of sig. corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
% 
% Posi_SigfCorrMean: mean of sig. positive corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
%
% Nega_SigfCorrMean: mean of sig. positive corr on each cell pair
% 1, initial cell;2-3, respond cell pair; 4, each corr
%
% NETWORK: cell pair matrix for each network
% 1:IntraIsland, 2:nonIsland, 3:InterIsland,

% cd(iwt_Dir{s,t});samplename=SampleName{3,s};
% caFr=10;
% cd(DIR)
dF= csvread(strcat(DIR, "\ST_PCI_noDup_dF.csv")); %
[NumFrame NumCell]=size(dF);
load(strcat(DIR,"\ST_dF_grid_aut_data.mat"), "lk");
%%
% trace 30s; sliding window 10s; 2s increment
INCRE=sWIN/5;
numINCRE=round((30-sWIN)/INCRE);

AA=nan(NumCell,NumCell);
for i=1:1:NumCell
    for ii=i+1:1:NumCell
        AA(i,ii)=1;
    end
end

SigfCorr_ratio=NaN(NumCell,NumCell,NumCell,9);% 1, initial cell;2-3, respond cell pair; 4, each corr
Posi_SigfCorrMean=NaN(NumCell,NumCell,NumCell);% 1, initial cell;2-3, respond cell pair
Nega_SigfCorrMean=NaN(NumCell,NumCell,NumCell);% 1, initial cell;2-3, respond cell pair
RR=[-0.8:0.2:0.8];
for c=1:1:NumCell

    % select Top 50 peak
    peak_dF=dF(round((lk{c})*caFr));
    [~, ID]=sort(peak_dF,"descend");
    if length(ID)>=50
        tmp_Num=50;
    else
       tmp_Num=length(ID);
    end
    rho=[];pval=[];
    k=1;
    for tmp_ID=1:1:tmp_Num; i=ID(tmp_ID);       
        for l=1:numINCRE
            tmp_FrST=round(lk{c}(i)*caFr+(l-1)*INCRE*caFr);
            tmp_FrED=round(tmp_FrST+sWIN*caFr);
            if tmp_FrED<NumFrame
                tmp_dF=dF(tmp_FrST:tmp_FrED,:);
                [r, p]=corr(tmp_dF);
                rho(:,:,k)=r.*AA;
                pval(:,:,k)=p.*AA;
                k=k+1;
            end
        end
    end

    for rr=1:1:9
        if RR(rr)<0
        tmp=pval<0.05&rho<RR(rr);
        SigfCorr_ratio(c,:,:,rr)=mean(tmp.*AA,3,"omitnan");% ratio of sig. corr on each cell pair
        elseif RR(rr)>0
        tmp=pval<0.05&rho>RR(rr);
        SigfCorr_ratio(c,:,:,rr)=mean(tmp.*AA,3,"omitnan");% ratio of sig. corr on each cell pair
        end
    end

    % mean of sig. positive corr on each cell pair
    tmp=nan(NumCell,NumCell,k-1);
    tmp(find(pval<0.05&rho>0))=1;
    Posi_SigfCorrMean(c,:,:)=mean((rho.*tmp),3,'omitnan');
    % mean of sig. negative corr on each cell pair
    tmp=nan(NumCell,NumCell,k-1);
    tmp(find(pval<0.05&rho<0))=1;
    Nega_SigfCorrMean(c,:,:)=mean((rho.*tmp),3,'omitnan');
end




end
