function [Minmax_Mean_ErrDist, Minmax_Mean_zCORRMAP, Minmean_Mean_ErrDist, Minmean_Mean_zCORRMAP]=Fn_CorrDecode_NormCellNum_230314(DIR)
close all;
% Mod = wt_Mod{s,t};
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')
% dF = csvread(strcat(DIR,"\ST_PCI_noDup_dF.csv")); %
% CellPos = csvread(strcat(DIR,"\ST_PCI_noDup_CellPos.csv")); 

x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))

CDir=pwd;



num_AnaCellset=50; % number of one cell set
num_SHUFFLE=20;% number of one cell set
num_SHUFFLE2=5;
numCells=size(GSrate_map,1);

Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 
vlim=2;
moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)


numFrame = size(Trk,1);
num_mFrame=size(moveFrame,1);
%%
if numCells<num_AnaCellset
    num_SHUFFLE=1
    num_SHUFFLE2=1
    num_AnaCellset=numCells
end

%% Frame set (7000 frames each)
num_RAND=ceil(num_mFrame/7000);
tmp_merge=ceil((num_RAND*7000-num_mFrame)/(num_RAND-1));
RAND=[];
for i=1:num_RAND
    RAND(i)=1+(6999-tmp_merge)*(i-1);
end


%%
for c=1:numCells
    TotalRatemap(c,:,:)=GSrate_map{c};
end

for SHUFFLE2=1:num_SHUFFLE2
    
    temp_cellset=nan(num_SHUFFLE,num_AnaCellset);
    for SHUFFLE = 1:num_SHUFFLE
        temp_cellset(SHUFFLE,:)=randsample(1:1:numCells,num_AnaCellset);
    end
    
    
    
    
    
    
    
    
    CORRMAP=nan(num_SHUFFLE,num_RAND,50,50,50,50);
    % tmpGSrate_map=nan(50,50);
    for ff=1:num_RAND
    
        mFRAME=RAND(ff):(RAND(ff)+7000) ;
        [tmpGSrate_map]=Fn_rate_map_forCorr_Mframe(DIR,mFRAME);
        for SHUFFLE = 1:num_SHUFFLE
            CELLSET=temp_cellset(SHUFFLE,:);
           
        
            for BaseY=1:50
                for y=1:50
                    tmp=corr(tmpGSrate_map(CELLSET,:,BaseY), ...
                        TotalRatemap(CELLSET,:,y), 'rows','complete');
                    for x=1:50
                        CORRMAP(SHUFFLE,ff,:,BaseY,x,y)=tmp(:,x);
                    end
                
                end
            end
        end
    end
    disp(strcat(num2str(SHUFFLE2),'Cormap end'))
    
    %{
    for ff=1:num_RAND
        for BaseX=1:10:50
            for BaseY=1:10:50 

                tmp(1:50,1:50)=CORRMAP(ff,BaseX,BaseY,:,:);
                imshow(tmp,'InitialMagnification',5000);
                text(1,2,strcat(samplename,", X=",num2str(BaseX),", Y=",num2str(BaseY)),"Color",'y')
                text(BaseY-0.5,BaseX+0.5,"*","Color","w","FontSize",17)
                colormap(jet)
                exportgraphics(gcf,strcat(CDir, "\Corr\each\" ,samplename," X",num2str(BaseX)," Y",num2str(BaseY)," f",num2str(ff),".jpg"))
                clf
        
            end
        end
    end
    %}
    %%
    
    ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    zCORRMAP=nan(num_SHUFFLE, num_RAND,50,50);
    tmp=[];
    
    for BaseX=1:50
        for BaseY=1:50 
            for ff=1:num_RAND
                for SHUFFLE = 1:num_SHUFFLE
                tmp(1:50,1:50)=CORRMAP(SHUFFLE,ff,BaseX,BaseY,:,:);
    
                    Maxpeak=max(tmp,[],'all');
                    [px py]=find(tmp==Maxpeak);
                    if isempty(px) ==0
                        ErrDist(SHUFFLE,ff,BaseX,BaseY)=sqrt((BaseX-px(1))^2 +(BaseY-py(1))^2)*2;
            
                        MEAN=mean( tmp,'all','omitnan');
                        STD=std(tmp,[],'all','omitnan');
                        POINT=tmp(BaseX,BaseY);
                        zCORRMAP(SHUFFLE, ff,BaseX,BaseY)=(POINT-MEAN)/STD;
        
                    end
                end
        
            end
            
        end
    end

    %{
    for ff=1:num_RAND
        for BaseX=1:10:50
            for BaseY=1:10:50 
if  ErrDist(SHUFFLE,ff,BaseX,BaseY) > 50

                tmp(1:50,1:50)=CORRMAP(SHUFFLE,ff,BaseX,BaseY,:,:);
                imshow(tmp,'InitialMagnification',5000);
                text(1,2,strcat(samplename,", X=",num2str(BaseX),", Y=",num2str(BaseY)),"Color",'y')
                text(BaseY-0.5,BaseX+0.5,"*","Color","w","FontSize",17)
                colormap(jet)
                exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\ex\Erro50 " ,samplename," X",num2str(BaseX)," Y",num2str(BaseY)," f",num2str(ff),".jpg"))
                clf
break
end
        
            end
        end
    end
    %}
    
    
    tmp=max(ErrDist,1,'omitnan');
    Mean_ErrDist(1:num_SHUFFLE,1:50,1:50)=mean(ErrDist,2,'omitnan');
    Mean_zCORRMAP(1:num_SHUFFLE,1:50,1:50)=mean(zCORRMAP,2,'omitnan');

    % minumum largest error
    tmp_ErrD=nan(num_SHUFFLE,1);
    for ss=1:num_SHUFFLE
        tmp_ErrD(ss)=max(Mean_ErrDist(ss,:,:),[],'all');
    end
    tmpCELLSET_id=find(tmp_ErrD==min(tmp_ErrD));
    
    tmp1Min_Mean_ErrDist(SHUFFLE2,1:50,1:50)=Mean_ErrDist(tmpCELLSET_id(1),:,:);
    tmp1Max_Mean_zCORRMAP(SHUFFLE2,1:50,1:50)=Mean_zCORRMAP(tmpCELLSET_id(1),:,:);
    
    % minumum mean error
    tmp_ErrD=nan(num_SHUFFLE,1);
    for ss=1:num_SHUFFLE
        tmp_ErrD(ss)=mean(Mean_ErrDist(ss,:,:),'all');
    end
    tmpCELLSET_id=find(tmp_ErrD==min(tmp_ErrD));
    
    tmp2Min_Mean_ErrDist(SHUFFLE2,1:50,1:50)=Mean_ErrDist(tmpCELLSET_id(1),:,:);
    tmp2Max_Mean_zCORRMAP(SHUFFLE2,1:50,1:50)=Mean_zCORRMAP(tmpCELLSET_id(1),:,:);




end

%%
% minumum largest error
tmp_ErrD=nan(num_SHUFFLE2,1);
for ss=1:num_SHUFFLE2
    tmp_ErrD(ss)=max(tmp1Min_Mean_ErrDist(ss,:,:),[],'all');
end
tmpCELLSET_id=find(tmp_ErrD==min(tmp_ErrD));

Minmax_Mean_ErrDist(1:50,1:50)=tmp1Min_Mean_ErrDist(tmpCELLSET_id(1),:,:);
Minmax_Mean_zCORRMAP(1:50,1:50)=tmp1Max_Mean_zCORRMAP(tmpCELLSET_id(1),:,:);

% minumum mean error
tmp_ErrD=nan(num_SHUFFLE2,1);
for ss=1:num_SHUFFLE2
    tmp_ErrD(ss)=mean(tmp2Min_Mean_ErrDist(ss,:,:),'all');
end
tmpCELLSET_id=find(tmp_ErrD==min(tmp_ErrD));

Minmean_Mean_ErrDist(1:50,1:50)=tmp1Min_Mean_ErrDist(tmpCELLSET_id(1),:,:);
Minmean_Mean_zCORRMAP(1:50,1:50)=tmp1Max_Mean_zCORRMAP(tmpCELLSET_id(1),:,:);


%%
imagesc(Minmax_Mean_ErrDist);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Mean Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest large error from 50 cell/100 shuffle")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\MeanError\Minmax_MeanError " ,samplename,".jpg"))
        clf
    

imagesc(Minmax_Mean_zCORRMAP);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Z Score of correlation",'FontSize',11);
title({samplename;strcat("Smallest large error from 50 cell/100 shuffle")},...
    'FontSize',14);

exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\MeanzCorr\Minmax_zCorr " ,samplename,".jpg"))
            clf

%%
imagesc(Minmean_Mean_ErrDist);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Mean Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest mean error from 50 cell/100 shuffle")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\MeanError\Minmean_MeanError " ,samplename,".jpg"))
        clf
    

imagesc(Minmax_Mean_zCORRMAP);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Z Score of correlation",'FontSize',11);
title({samplename;strcat("Smallest mean error from 50 cell/100 shuffle")},...
    'FontSize',14);

exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\MeanzCorr\Minmean_zCorr " ,samplename,".jpg"))
            clf




end



