function [MEAN, TOTAL, MAX]=Fn_CorrDecode_NC_zRatemap_230804(DIR)
close all;
% Mod = wt_Mod{s,t};
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map','GSoccup_m')
% dF = csvread(strcat(DIR,"\ST_PCI_noDup_dF.csv")); %
% CellPos = csvread(strcat(DIR,"\ST_PCI_noDup_CellPos.csv")); 

x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))

CDir=pwd;



num_AnaCellset=50; % number of one cell set
num_SHUFFLE=50;% number of one cell set
num_SHUFFLE2=20;
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


%% Z scored rate map
TotalRatemap=nan(numCells,50,50);
for c=1:numCells
    tmp=GSrate_map{c}; 
    tmp(GSoccup_m<10)= nan; % remove bin less than 10 frame
    TotalRatemap(c,:,:)=tmp;
    
end

Z_TotalRatemap=nan(numCells,50,50);
MEAN=mean(TotalRatemap,[2,3],"omitnan");
STD=std(TotalRatemap,[],[2 3],"omitnan");
Z_TotalRatemap=(TotalRatemap-MEAN)./STD;


%% make matrix for non specific case
x=   2:5:50;
y=   2:5:50;
k=1;
for ix=1:length(x)
    for iy=1:length(x)
        NONSPECIC(k,:)=[x(ix) y(iy)];
        k=k+1;

    end
end
%% calculate error distnce on each frame set, cell set, each spatial location(bin)

for SHUFFLE2=1:num_SHUFFLE2
    
    temp_cellset=nan(num_SHUFFLE,num_AnaCellset);
    for SHUFFLE = 1:num_SHUFFLE
        temp_cellset(SHUFFLE,:)=randsample(1:1:numCells,num_AnaCellset);
    end
    

    CORRMAP=nan(num_SHUFFLE,num_RAND,50,50,50,50);
    % tmpGSrate_map=nan(50,50);
    for ff=1:num_RAND
    
        mFRAME=RAND(ff):(RAND(ff)+7000) ;
        [Z_tmpGSrate_map]=Fn_zRatemap_forCorr_Mframe(DIR,mFRAME);
% xxxtmp=[];
% xxxtmp(1:50,1:50)=Z_tmpGSrate_map(2,:,:);
% imshow(xxxtmp,'InitialMagnification',5000);
        for SHUFFLE = 1:num_SHUFFLE
            CELLSET=temp_cellset(SHUFFLE,:);
           
        
            for BaseY=1:50
                for y=1:50
                    tmp=corr(Z_tmpGSrate_map(CELLSET,:,BaseY), ...
                        Z_TotalRatemap(CELLSET,:,y), 'rows','pairwise');
                    % if isnan(tmp)==1
                    %     BaseY
                    %     y
                    % 
                    % end

                    for x=1:50
                        CORRMAP(SHUFFLE,ff,:,BaseY,x,y)=tmp(:,x);
                    end
                
                end
            end
        
% xxxtmp=[];
% xxxtmp(1:50,1:50)=CORRMAP(SHUFFLE,ff,20,20,:,:);
% imshow(xxxtmp,'InitialMagnification',5000)
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
    
    zCORRMAP=nan(num_SHUFFLE, num_RAND,50,50,50,50);
    mean_ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    total_ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    max_ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    

    for BaseX=1:50
        for BaseY=1:50 
            for ff=1:num_RAND
                for SHUFFLE = 1:num_SHUFFLE
                    tmp=[];
                    tmp(1:50,1:50)=CORRMAP(SHUFFLE,ff,BaseX,BaseY,:,:);
                    MEAN=mean( tmp,'all','omitnan');
                    STD=std(tmp,[],'all','omitnan');
              % POINT=tmp(BaseX,BaseY);
                    zCORRMAP(SHUFFLE, ff,BaseX,BaseY,:,:)=(tmp-MEAN)/STD;

                    z_tmp=[];
                    z_tmp(1:50,1:50)=zCORRMAP(SHUFFLE,ff,BaseX,BaseY,:,:);
                    % imshow(z_tmp,'InitialMagnification',5000);
                % find peak in similarity map (z score>=2)
                    I2bi=zeros(50,50);%made binary image
                    I2bi(z_tmp>=2)=1 ;%made binary image
                    % imshow(I2bi,'InitialMagnification',5000);
                
                % Mark position of centroid 
                    bwl = bwlabel(I2bi); %label binary object
                    biarea = regionprops(bwl, 'Area', 'BoundingBox', 'Centroid'); 
                    tmp_err=nan(1);
                    if  isempty(vertcat(biarea.Area)) ==0     
                        area_array = vertcat(biarea.Area);
                        center_array = vertcat(biarea.Centroid); 
                        center_array(area_array<8,:)=[];; % cut small error peak (3x3)
                 % calculate error distance to peaks 
                 %%% CUTION! center_array =(y,x)
                        tmp_err=pdist2([BaseX BaseY], [center_array(:,2) center_array(:,1)]);% calculate error distance 
                
                    end
                     % if  isnan(tmp_err) ==1  ;% 
                     %     tmp_err=pdist2([BaseX BaseY], NONSPECIC);
                     % elseif isempty(tmp_err) ==1  ;% no specific point
                     if isempty(tmp_err) ==1  ;% no specific peak in zCORRMAP
                          tmp_err=pdist2([BaseX BaseY], NONSPECIC);
                     end

                    mean_ErrDist(SHUFFLE,ff,BaseX,BaseY)=mean(tmp_err,"omitnan");% MEAN error distance
                    total_ErrDist(SHUFFLE,ff,BaseX,BaseY)=sum(tmp_err,"omitnan");% 
                    max_ErrDist(SHUFFLE,ff,BaseX,BaseY)=max(tmp_err,[],"omitnan");% max error distance
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
    
    


% average result from rand frame set
    Mean_ErrDist_meanframeset(1:num_SHUFFLE,1:50,1:50)=mean(mean_ErrDist,2,'omitnan'); % 
    Total_ErrDist_meanframeset(1:num_SHUFFLE,1:50,1:50)=mean(total_ErrDist,2,'omitnan'); % average result from rand frame set
    Max_ErrDist_meanframeset(1:num_SHUFFLE,1:50,1:50)=mean(max_ErrDist,2,'omitnan'); % average result from rand frame set
    Mean_zCORRMAP_meanframeset(1:num_SHUFFLE,1:50,1:50,1:50,1:50)=mean(zCORRMAP,2,'omitnan'); % average result from rand frame set

    % tmp=[];
    % tmp(1:50,1:50)=mean_ErrDist(1,1,:,:);
    % imagesc(tmp);

% chose cellset with minumum total (mean) error on each spatial bin
    tmp_Mean_ErrDist_meanframeset_bestcellset(SHUFFLE2,1:50,1:50)=min(Mean_ErrDist_meanframeset,[],1);
    tmp_Total_ErrDist_meanframeset_bestcellset(SHUFFLE2,1:50,1:50)=min(Total_ErrDist_meanframeset,[],1);
    tmp_Max_ErrDist_meanframeset_bestcellset(SHUFFLE2,1:50,1:50)=min(Max_ErrDist_meanframeset,[],1);

end


% chose cellset with minumum total (mean) error on each spatial bin
    Mean_ErrDist_meanframeset_bestcellset(1:50,1:50)=min(tmp_Mean_ErrDist_meanframeset_bestcellset,[],1);
    Total_ErrDist_meanframeset_bestcellset(1:50,1:50)=min(tmp_Total_ErrDist_meanframeset_bestcellset,[],1);
    Max_ErrDist_meanframeset_bestcellset(1:50,1:50)=min(tmp_Max_ErrDist_meanframeset_bestcellset,[],1);

    MEAN=Mean_ErrDist_meanframeset_bestcellset;
    TOTAL=Total_ErrDist_meanframeset_bestcellset;
    MAX=Max_ErrDist_meanframeset_bestcellset;

% 
% Minmean_Mean_ErrDist(1:50,1:50)=tmp1Min_Mean_ErrDist(tmpCELLSET_id(1),:,:);
% Minmean_Mean_zCORRMAP(1:50,1:50)=tmp1Max_Mean_zCORRMAP(tmpCELLSET_id(1),:,:);


%% visualise
imagesc(Mean_ErrDist_meanframeset_bestcellset);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Mean Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest mean error from 50 cell/",num2str(num_SHUFFLE*num_SHUFFLE2)," shuffle")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\zMeanError\zMeanError averagedFrameSet bestCellSet" ,samplename,".jpg"))
        clf
    

imagesc(Total_ErrDist_meanframeset_bestcellset);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Total Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest total error from 50 cell/",num2str(num_SHUFFLE*num_SHUFFLE2)," shuffle")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\zTotalError\zTotalError averagedFrameSet bestCellSet" ,samplename,".jpg"))
        clf
    
imagesc(Max_ErrDist_meanframeset_bestcellset);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Max Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest Max error from 50 cell/",num2str(num_SHUFFLE*num_SHUFFLE2)," shuffle")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\zMaxError\zTotalError averagedFrameSet bestCellSet" ,samplename,".jpg"))
        clf
    



end



