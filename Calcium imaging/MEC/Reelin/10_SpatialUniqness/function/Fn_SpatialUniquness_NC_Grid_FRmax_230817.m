function [MEAN, TOTAL, MAX]=Fn_SpatialUniquness_NC_Grid_FRmax_230817(DIR,caFr)
close all;
% Mod = wt_Mod{s,t};
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map','GSoccup_m','Grid_Cells')
% dF = csvread(strcat(DIR,"\ST_PCI_noDup_dF.csv")); %
% CellPos = csvread(strcat(DIR,"\ST_PCI_noDup_CellPos.csv")); 

xxx=split(DIR,"\");
samplename=strcat(char(strrep(xxx(end-1),'_','-')), " ",char(strrep(strrep(xxx(end),'_',' '),'OF','')))
xxx=[];

CDir=pwd;

% numCells=size(GSrate_map,1);
num_SelectCell=10; % number of select cell set
AnaCellset=Grid_Cells;
num_SelectFrame=caFr*1000;
num_SHUFFLE=1;


Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 
vlim=2;
moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)


numFrame = size(Trk,1);
num_mFrame=size(moveFrame,1);
%% Make cell set to analyise

if num_SelectCell > AnaCellset
    num_SelectCell = AnaCellset
end


%% Frame set (15000 frames each)
num_RAND=ceil(num_mFrame/num_SelectFrame);
tmp_merge=ceil((num_RAND*num_SelectFrame-num_mFrame)/(num_RAND-1));
RAND=[];
for i=1:num_RAND
    RAND(i)=1+(num_SelectFrame-1-tmp_merge)*(i-1);
end




% %% make matrix for non specific case
% x=   2:5:50;
% y=   2:5:50;
% k=1;
% for ix=1:length(x)
%     for iy=1:length(x)
%         NONSPECIC(k,:)=[x(ix) y(iy)];
%         k=k+1;
% 
%     end
% end

%
%% calculate error distnce on each frame set, cell set, each spatial location(bin)

         CORRMAP=nan(num_SHUFFLE,num_RAND,50,50,50,50);
         CORRMAP_Pval=nan(num_SHUFFLE,num_RAND,50,50,50,50);
tic
    for ff=1:num_RAND
        mFRAME=RAND(ff):(RAND(ff)+num_SelectFrame) ;
        tmpGSrate_map=nan(50,50);
        % [Z_tmpGSrate_map]=Fn_zRatemap_forCorr_Mframe(DIR,mFRAME);

        [tmpGSrate_map]=Fn_rate_map_forCorr_Mframe(DIR,mFRAME);
      
        parfor (BaseX=1:50,2)
        % for BaseX=1:50
            for BaseY=1:50              
                [CORRMAP(:,ff,BaseX,BaseY,:,:) CORRMAP_Pval(:,ff,BaseX,BaseY,:,:)]=Fn_maxFR_CORRMAP_Grid(num_SHUFFLE,BaseX,BaseY,num_SelectCell,AnaCellset,tmpGSrate_map);
            end
        end   
    end
  toc
%tmp=reshape(CORRMAP(1,1,1,1,:,:),[],50);
%imshow(tmp,"InitialMagnification",5000)

    %%
    
    zCORRMAP=nan(num_SHUFFLE, num_RAND,50,50,50,50);
    mean_ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    total_ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    max_ErrDist=nan(num_SHUFFLE, num_RAND,50,50);
    

    parfor (BaseX=1:50,2)
    % for BaseX=1:50
        for BaseY=1:50 
            for ff=1:num_RAND
                for SHUFFLE = 1:num_SHUFFLE
                    tmp=[];tmp2=[];
                    tmp(1:50,1:50)=CORRMAP(SHUFFLE,ff,BaseX,BaseY,:,:);
                    tmp2(1:50,1:50)=CORRMAP_Pval(SHUFFLE,ff,BaseX,BaseY,:,:);
                    
                    %  imshow(tmp*-1,'InitialMagnification',5000);
%{
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
                imshow(I2bi,'InitialMagnification',5000);
%}

          % find peak in similarity map (r>=0.2&pval<0.01)
                I2bi=zeros(50,50);%made binary image
                    I2bi(tmp>=0.2&tmp2<0.01)=1 ;%made binary image
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
                     if isempty(tmp_err) ==1  ;
                        area_array = vertcat(biarea.Area);
                        center_array = vertcat(biarea.Centroid); 
                        tmp_err=pdist2([BaseX BaseY], [center_array(:,2) center_array(:,1)]);% calculate error distance 
                        tmp_err=min(tmp_err);

                     end

                    mean_ErrDist(SHUFFLE,ff,BaseX,BaseY)=mean(tmp_err,"omitnan");% MEAN error distance
                    total_ErrDist(SHUFFLE,ff,BaseX,BaseY)=sum(tmp_err,"omitnan");% 
                    max_ErrDist(SHUFFLE,ff,BaseX,BaseY)=max(tmp_err,[],"omitnan");% max error distance
                end
        
            end
            
        end
    end
% toc
    

% average result from rand frame set
    Mean_ErrDist_meanframeset(1:num_SHUFFLE,1:50,1:50)=mean(mean_ErrDist,2,'omitnan'); % 
    Total_ErrDist_meanframeset(1:num_SHUFFLE,1:50,1:50)=mean(total_ErrDist,2,'omitnan'); % average result from rand frame set
    Max_ErrDist_meanframeset(1:num_SHUFFLE,1:50,1:50)=mean(max_ErrDist,2,'omitnan'); % average result from rand frame set
    % Mean_zCORRMAP_meanframeset(1:num_SHUFFLE,1:50,1:50,1:50,1:50)=mean(zCORRMAP,2,'omitnan'); % average result from rand frame set






% chose cellset with minumum total (mean) error on each spatial bin 
    Mean_ErrDist_meanframeset_bestcellset(1:50,1:50)=min(Mean_ErrDist_meanframeset,[],1);
    Total_ErrDist_meanframeset_bestcellset(1:50,1:50)=min(Total_ErrDist_meanframeset,[],1);
    Max_ErrDist_meanframeset_bestcellset(1:50,1:50)=min(Max_ErrDist_meanframeset,[],1);

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
strcat("Mean error from Top ",num2str(num_SelectCell)," FR Grid cell")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\MeanError_GridmaxFR\Grid MeanError averagedFrameSet maxFR " ,samplename,".jpg"))
        clf
    
% tmp=Total_ErrDist_meanframeset_bestcellset;
% tmp(find(tmp>20))=nan;
% imagesc(tmp)


imagesc(Total_ErrDist_meanframeset_bestcellset);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Total Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Total error from Top ",num2str(num_SelectCell)," FR Grid cell")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\TotalError_GridmaxFR\Grid TotalError averagedFrameSet maxFR " ,samplename,".jpg"))
        clf
    
imagesc(Max_ErrDist_meanframeset_bestcellset);
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet);
cb=colorbar;

ylabel(cb,"Max Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Max error from Top ",num2str(num_SelectCell)," FR Grid cell")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\MaxError_GridmaxFR\Grid MaxError averagedFrameSet maxFR " ,samplename,".jpg"))
        clf
    



end



