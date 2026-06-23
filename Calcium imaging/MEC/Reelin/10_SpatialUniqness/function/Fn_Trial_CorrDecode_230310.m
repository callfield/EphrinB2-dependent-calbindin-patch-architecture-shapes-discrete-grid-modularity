function Fn_Trial_CorrDecode_230310(DIR)


close all;
% Mod = wt_Mod{s,t};
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')
% dF = csvread(strcat(DIR,"\ST_PCI_noDup_dF.csv")); %
% CellPos = csvread(strcat(DIR,"\ST_PCI_noDup_CellPos.csv")); 

x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))

CDir=pwd;


%%

for c=1:size(GSrate_map,1)
    TotalRatemap(c,:,:)=GSrate_map{c};
end


Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 
numFrame = size(Trk,1);

ff=1
    FRAME=1:numFrame;
    [tmpGSrate_map]=Fn_rate_map_forCorr(DIR,FRAME);

    for BaseY=1:50
        for y=1:50
            tmp=corr(tmpGSrate_map(:,:,BaseY), TotalRatemap(:,:,y), 'rows','complete');
            for x=1:50
                CORRMAP(ff,:,BaseY,x,y)=tmp(:,x);
            end
        
        end
    end

    for BaseX=1:10:50
        for BaseY=1:10:50 
            tmp(1:50,1:50)=CORRMAP(ff,BaseX,BaseY,:,:);
            imshow(tmp,'InitialMagnification',5000);
            text(1,2,strcat(samplename,", X=",num2str(BaseX),", Y=",num2str(BaseY)),"Color",'y')
            text(BaseY-0.5,BaseX+0.5,"*","Color","w","FontSize",17)
            colormap(jet)
            exportgraphics(gcf,strcat(CDir, "\Corr\each\" ,samplename," X",num2str(BaseX)," Y",num2str(BaseY),".jpg"))
            clf
    
        end
    end


end