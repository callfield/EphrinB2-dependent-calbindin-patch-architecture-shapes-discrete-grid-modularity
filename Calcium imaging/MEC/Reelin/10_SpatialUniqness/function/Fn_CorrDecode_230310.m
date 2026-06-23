function [Mean_ErrDist, Mean_zCORRMAP,ErrDist,zCORRMAP]=Fn_CorrDecode_230310(DIR)
close all;
% Mod = wt_Mod{s,t};
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')
% dF = csvread(strcat(DIR,"\ST_PCI_noDup_dF.csv")); %
% CellPos = csvread(strcat(DIR,"\ST_PCI_noDup_CellPos.csv")); 

x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))

CDir=pwd;



% 



%%

for c=1:size(GSrate_map,1)
    TotalRatemap(c,:,:)=GSrate_map{c};
end


Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 
numFrame = size(Trk,1);

RAND=10000:10000:numFrame;
RAND(length(RAND)+1)=numFrame;
CORRMAP=nan(length(RAND),50,50,50,50);
tmpGSrate_map=nan(50,50);
for ff=1:length(RAND)

    FRAME=(RAND(ff)-9999):RAND(ff);
    [tmpGSrate_map]=Fn_rate_map_forCorr(DIR,FRAME);

    for BaseY=1:50
        for y=1:50
            tmp=corr(tmpGSrate_map(:,:,BaseY), TotalRatemap(:,:,y), 'rows','complete');
            for x=1:50
                CORRMAP(ff,:,BaseY,x,y)=tmp(:,x);
            end
        
        end
    end
end

%{
for ff=1:length(RAND)
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

ErrDist=nan(length(RAND),50,50);
zCORRMAP=nan(length(RAND),50,50);
tmp=[];

for BaseX=1:50
    for BaseY=1:50 
        for ff=1:length(RAND)
            tmp(1:50,1:50)=CORRMAP(ff,BaseX,BaseY,:,:);

            Maxpeak=max(tmp,[],'all');
            [px py]=find(tmp==Maxpeak);
            if isempty(px) ==0
                ErrDist(ff,BaseX,BaseY)=sqrt((BaseX-px)^2 +(BaseY-py)^2)*2;
    
                MEAN=mean( tmp,'all','omitnan');
                STD=std(tmp,[],'all','omitnan');
                POINT=tmp(BaseX,BaseY);
                zCORRMAP(ff,BaseX,BaseY)=(POINT-MEAN)/STD;

            end
    
        end
        
    end
end

Mean_ErrDist(1:50,1:50)=mean(ErrDist,1,'omitnan');
Mean_zCORRMAP(1:50,1:50)=mean(zCORRMAP,1,'omitnan');


%%
imagesc(Mean_ErrDist);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Mean Error Distance (cm)",'FontSize',11);
title(strcat( samplename ),'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\Corr\MeanError\MeanError " ,samplename,".jpg"))
        clf
    

imagesc(Mean_zCORRMAP);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Z Score of correlation",'FontSize',11);
title(strcat( samplename ),'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\Corr\MeanzCorr\zCorr " ,samplename,".jpg"))
            clf




end



