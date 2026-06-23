function [ Min_ErrDist]=Fn_CorrDecode_NormCellNum_GridMod_230322(DIR,All_CELLSET,NAME)
close all;

% All_CELLSET=Grid_Cells(wt_Mod{s,t}{1});
load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'GSrate_map')


x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_','-')), " ",char(strrep(strrep(x(end),'_',' '),'OF','')))

CDir=pwd;



num_AnaCellset=10; % number of one cell set
num_AnaFrameset=15000; % number of one frame set
num_SHUFFLE=20;% number of make cell set (Total shuffle = num_SHUFFLE x num_SHUFFLE2)
% num_SHUFFLE2=50;
numCells=size(All_CELLSET,2);



Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 
vlim=2;
moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)


numFrame = size(Trk,1);
num_mFrame=size(moveFrame,1);
%%
if numCells<=num_AnaCellset
    num_SHUFFLE=1
    num_SHUFFLE2=1
    num_AnaCellset=numCells
end

if num_mFrame<num_AnaFrameset
    num_AnaFrameset=num_mFrame
end

%% Frame set ("7000" frames each)
num_FrameSet=ceil(num_mFrame/num_AnaFrameset);
tmp_merge=ceil((num_FrameSet*num_AnaFrameset-num_mFrame)/(num_FrameSet-1));
RAND=[];
for i=1:num_FrameSet
    RAND(i)=1+(num_AnaFrameset-tmp_merge-1)*(i-1);
end

if num_FrameSet==1
    RAND=1;
end



NC_CELLSET=nchoosek(All_CELLSET,num_AnaCellset); % all possible cell set

num_SHUFFLE2=ceil(size(NC_CELLSET,1)/num_SHUFFLE);
% if num_SHUFFLE2 is too big, change to random sampling
if num_SHUFFLE2>50
    num_SHUFFLE2=50
    tmp_cell=randsample(1:1:size(NC_CELLSET,1),num_SHUFFLE2*num_SHUFFLE);
    NC_CELLSET=NC_CELLSET(tmp_cell,:);
else
    num_SHUFFLE2
end




for SHUFFLE2=1:num_SHUFFLE2
    
    st=1+num_SHUFFLE*(SHUFFLE2-1);
    ed=num_SHUFFLE*(SHUFFLE2);
    if ed > size(NC_CELLSET,1)
        ed=size(NC_CELLSET,1);
        num_SHUFFLE=ed-st+1;
    end
    temp_cellset=NC_CELLSET(st:ed,:);
    
    CORRMAP=nan(num_SHUFFLE,num_FrameSet,50,50,50,50);
    % tmpGSrate_map=nan(50,50);
    for ff=1:num_FrameSet
    
        mFRAME=RAND(ff):(RAND(ff)+num_AnaFrameset-1) ;
        [tmpGSrate_map]=Fn_rate_map_forCorr_Mframe(DIR,mFRAME);
        for SHUFFLE = 1:num_SHUFFLE
            CELLSET=temp_cellset(SHUFFLE,:);
           
            for BaseY=1:50
                for y=1:50

                    tmp=corr(tmpGSrate_map(CELLSET,:,BaseY), ...
                        tmpGSrate_map(CELLSET,:,y), 'rows','complete');
                    for x=1:50
                        CORRMAP(SHUFFLE,ff,:,BaseY,x,y)=tmp(:,x);
                    end
                
                end
            end

        end
    end


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
    

    Mean_CORRMAP=mean(CORRMAP,2,"omitnan"); % average difference on frame set
    CORRMAP=[];

    ErrDist=nan(num_SHUFFLE,50,50);
    tmp=[];
    
    for BaseX=1:50
        for BaseY=1:50 
            for SHUFFLE = 1:num_SHUFFLE

                tmp=reshape(Mean_CORRMAP(SHUFFLE,1,BaseX,BaseY,:,:),50,50);
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
                    
                    DIST=sqrt(sum(([BaseY, BaseX]-XY_tmp).^2,2))*2; %  caution! XY is opposite
                    ErrDist(SHUFFLE,BaseX,BaseY)=mean(DIST);
         
                end
            end
        end
    end
 
    % chose cell set with minumum error dist on each bin
    tmpMin_ErrDist(SHUFFLE2,1:50,1:50)=min(ErrDist,[],1,"omitnan");

end

% minumum error
Min_ErrDist=reshape(min(tmpMin_ErrDist,[],1,"omitnan"),50,50);


%%
imagesc(Min_ErrDist,'AlphaData',~isnan(Min_ErrDist));% nan to transparent
set(gca,'color', [1 1 1]); % set backgound color
xticks([]); yticks([]);
daspect([100 100 100])
colormap(jet); 
cb=colorbar;

ylabel(cb,"Mean Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest error from ",num2str(num_AnaCellset), ...
     " grid cell with same module")},...
    'FontSize',14);



ylabel(cb,"Error Distance (cm)",'FontSize',11);
title({samplename;...
strcat("Smallest error from ",num2str(num_AnaCellset), ...
     " grid cell with same module")},...
    'FontSize',14);
exportgraphics(gcf,strcat(CDir, "\Corr_NormCellNum\GridModule\MinError\",NAME," Minmean_MeanError " ,samplename,".jpg"))
        clf
    
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



