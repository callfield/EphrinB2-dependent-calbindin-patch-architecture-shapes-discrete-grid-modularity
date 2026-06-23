function [ Min_meanErrDist Min_maxErrDist]=Fn_CorrDecode_NormCellNum_GridMod_230324(DIR,All_CELLSET,R2ID,NAME)
close all;
x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))
x=[];
% load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'Grid_Cells');
% All_CELLSET=Grid_Cells;
% All_CELLSET=Grid_Cells(wt_Mod{s,t}{1});

CDir=pwd;

num_AnaCell=10; % number of one cell set
num_SelectCell=12; % select Top12 contributed cells on each bin.
num_AnaFrameset=15000; % number of one frame set
num_SHUFFLE=33;
% num_SHUFFLE2=2;
numCells= max(size(All_CELLSET));

Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 
vlim=2;
moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
numFrame = size(Trk,1);
num_mFrame=size(moveFrame,1);



if numCells > num_AnaCell


%% Make frame set to analyise
if num_AnaFrameset>num_mFrame
    num_AnaFrameset=num_mFrame;
end


num_FrameSet=ceil(num_mFrame/num_AnaFrameset);
tmp_merge=ceil((num_FrameSet*num_AnaFrameset-num_mFrame)/(num_FrameSet-1));
Frame_st=[];
for i=1:num_FrameSet
    Frame_st(i)=1+(num_AnaFrameset-tmp_merge-1)*(i-1);
end

if num_FrameSet==1
    Frame_st=1;
end



%% Make cell set to analyise

if num_SelectCell > numCells
    num_SelectCell = numCells
end

NC_CELLSET=cell(50,50);
for BaseX=1:50
    for BaseY=1:50
        [C,ia,ib] = intersect(R2ID{BaseX,BaseY},All_CELLSET);
        % ia: rank of R^2 (1 is highest)
        [sort_ia sort_iaID] = sort(ia,'ascend');
        sorted_All_CELLSET=All_CELLSET(ib(sort_iaID)); % All_CELLSET sorted with rank of R^2
        
        tmp_CELLSET=sorted_All_CELLSET(1:num_SelectCell);% select Top12 contributed cells on each bin.
        NC_CELLSET{BaseX,BaseY}=nchoosek(tmp_CELLSET,num_AnaCell); % all possible cell set

    end
end


%% Number of shaffle
num_SHUFFLE2=ceil(size(NC_CELLSET{BaseX,BaseY},1)/num_SHUFFLE);


%%
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')
tmpGSrate_map=nan(numCells,50,50);
for c=1:numCells
    tmpGSrate_map(c,:,:)=GSrate_map{c};
end

% z scored
TotalRatemap=nan(numCells,50,50);
    MEAN=mean(tmpGSrate_map,[2,3],"omitnan");
    STD=std(tmpGSrate_map,[],[2 3],"omitnan");
    TotalRatemap=(tmpGSrate_map-MEAN)./STD;


%%

tmpMin_meanErrDist=[];
tmpMin_maxErrDist=[];
for SHUFFLE2=1:num_SHUFFLE2
    CORRMAP=nan(num_SHUFFLE,num_FrameSet,50,50,50,50);
    
    if SHUFFLE2==num_SHUFFLE2
        tmp=mod(size(NC_CELLSET{BaseX,BaseY},1),num_SHUFFLE);
        if tmp>0
            num_SHUFFLE=mod(size(NC_CELLSET{BaseX,BaseY},1),num_SHUFFLE);
        end
    end
        
        
%     tic
    for ff=1:num_FrameSet
        mFRAME=Frame_st(ff):(Frame_st(ff)+num_AnaFrameset-1) ;
        [tmpGSrate_map]=Fn_rate_map_forCorr_Mframe(DIR,mFRAME);
        MEAN=mean(tmpGSrate_map,[2,3],"omitnan");
        STD=std(tmpGSrate_map,[],[2 3],"omitnan");
        tmp_zGSrate_map=(tmpGSrate_map-MEAN)./STD;
        for BaseX=1:50
%             tic
            for BaseY=1:50
                for SHUFFLE = 1:num_SHUFFLE
                   tmp_cellset=SHUFFLE + num_SHUFFLE*(SHUFFLE2-1);
                   CELLSET=NC_CELLSET{BaseX,BaseY}(tmp_cellset,:);
                   for x=1:50    
                       tmp_total=reshape( TotalRatemap(CELLSET,x,:),10, 50);
                       tmp=corr(tmp_zGSrate_map(CELLSET,BaseX,BaseY),...
                          tmp_total, 'rows','complete');
                       CORRMAP(SHUFFLE,ff,BaseX,BaseY,x,:)=tmp;
                   end
                end
            end
%             toc
        end
    end

% toc


%     disp(strcat(num2str(SHUFFLE2),' Cormap end'))
    %{
    for ff=1:num_FrameSet
        for BaseX=1:10:50
            for BaseY=1:10:50 

                tmp=reshape(CORRMAP(SHUFFLE,ff,BaseX,BaseY,:,:),50,50);
                imagesc(tmp)
                MEAN=mean( tmp,'all','omitnan');
                STD=std(tmp,[],'all','omitnan');
                
                tmp2=(tmp-MEAN)/STD;
                tmp3=zeros(50,50);
                tmp3(find(tmp2>2))=1;
                imagesc(tmp3)
                text(1,2,strcat(samplename,", X=",num2str(BaseX),", Y=",num2str(BaseY)),"Color",'y')
                text(BaseY-0.5,BaseX+0.5,"*","Color","w","FontSize",17)
                colormap(jet)
                exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\\GridModule\ex\" ,samplename," X",num2str(BaseX)," Y",num2str(BaseY)," f",num2str(ff),".jpg"))
                clf
        
            end
        end
    end
    %}
    

    ffMean_CORRMAP=mean(CORRMAP,2,"omitnan"); % average difference on frame set
%     CORRMAP=[];

    meanErrDist=nan(num_SHUFFLE,50,50);
    maxErrDist=nan(num_SHUFFLE,50,50);
    tmp=[];
    for BaseX=1:50
%         tic
        for BaseY=1:50 
            for SHUFFLE = 1:num_SHUFFLE

                tmp=reshape(ffMean_CORRMAP(SHUFFLE,1,BaseX,BaseY,:,:),50,50);
%                 imagesc(tmp)
                % find zscore (corr) > 2
                MEAN=mean( tmp,'all','omitnan');
                STD=std(tmp,[],'all','omitnan');                 
                tmp2=(tmp-MEAN)/STD;
                tmp3=zeros(50,50);
                 if isempty(find(tmp2>2)) ==0
                    tmp3(find(tmp2>2))=1;

                    %mark area, zscore(corr) > 2
                    bwl = bwlabel(tmp3); %label binary object
                    biarea = regionprops(bwl, 'Area', 'BoundingBox', 'Centroid'); %#ok<MRPBW>
                    area_array = vertcat(biarea.Area);
                    XY_array = vertcat(biarea.Centroid);
                    XY_tmp=XY_array(area_array>3,:); %only area more than 3 (4x4cm)
                    if isempty(XY_tmp)==0
                        tmp=sum(([BaseY, BaseX]-XY_tmp).^2,2);%  caution! XY is opposite
                        DIST=sqrt(tmp)*2; % bin=2cm

                        meanErrDist(SHUFFLE,BaseX,BaseY)=mean(DIST);% average of error on a cell set
                        maxErrDist(SHUFFLE,BaseX,BaseY)=max(DIST);% larget error on a cell set
                    end
                end
            end
        end
%         toc
    end
 
    % chose cell set with minumum (mean)error dist on each bin
    tmpMin_meanErrDist(SHUFFLE2,1:50,1:50)=min(meanErrDist,[],1,"omitnan");
    % chose cell set with minumum (max)error dist on each bin
    tmpMin_maxErrDist(SHUFFLE2,1:50,1:50)=min(maxErrDist,[],1,"omitnan");
end

Min_meanErrDist=[];
Min_maxErrDist=[];
% minumum error
Min_meanErrDist(1:50,1:50)=min(tmpMin_meanErrDist,[],1,"omitnan");
Min_maxErrDist(1:50,1:50)=min(tmpMin_maxErrDist,[],1,"omitnan");

%%
imagesc(Min_meanErrDist,'AlphaData',~isnan(Min_meanErrDist));% nan to transparent
set(gca,'color', [1 1 1]); % set backgound color
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet); 
cb=colorbar;

ylabel(cb,"Min average error distance (cm)",'FontSize',11);
title({samplename;...
strcat("Min average error (",num2str(num_AnaCell), ...
     " cell from ",NAME,")")},...
    'FontSize',14);

exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\GridModule\Min_meanError\",NAME," Min_meanError " ,samplename,".jpg"))
        clf
    
        
imagesc(Min_maxErrDist,'AlphaData',~isnan(Min_maxErrDist));% nan to transparent
set(gca,'color', [1 1 1]); % set backgound color
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet); 
cb=colorbar;

ylabel(cb,"Min largest error distance (cm)",'FontSize',11);
title({samplename;...
strcat("Min largest error(",num2str(num_AnaCell), ...
     " cell from ",NAME,")")},...
    'FontSize',14);

exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\GridModule\Min_maxError\",NAME," Min_maxError " ,samplename,".jpg"))
        clf
else
    Min_meanErrDist=nan(50,50);
    Min_maxErrDist=nan(50,50);
end
% 
% imagesc(Minmax_Mean_zCORRMAP);
% xticks([]); yticks([]);
% daspect([100 100 100])
% colormap(jet);
% cb=colorbar;
% 
% ylabel(cb,"Z Score of correlation",'FontSize',11);
% title({samplename; ...
%     strcat("Smallest large error from ",num2str(num_AnaCellset), ...
%     " cell/",num2str(num_SHUFFLE*num_SHUFFLE2)," shuffle")},...
%     'FontSize',14);
% 
% exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\GridModule\MeanzCorr\",NAME," Minmean_zCorr " ,samplename,".jpg"))
%             clf
% 

end



