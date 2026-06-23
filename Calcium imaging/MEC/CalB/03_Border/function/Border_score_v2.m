function [B_Score_v2]=Border_score_v2(CDIR,caFr)
close all;
load("ST_dF_grid_aut_data.mat")

x=split(pwd,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')))

s_dwSample=4;
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
STD=2; %Definition of peak
bin = 2 ; %Analysis each (cm)
win = 5 ; % moving window size
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);

Trk = csvread("ST_PCI_Ca_behav_track.csv"); % 
dF = csvread("ST_PCI_noDup_dF.csv"); %
[numFrames numCells] = size(dF);
CellPos = csvread("ST_PCI_noDup_CellPos.csv"); 
bg2=imresize(bg, s_dwSample, 'nearest');
vlim=2;


moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
move_dF = dF(moveFrame,:);
move_Trk = Trk(moveFrame,:);

%% Define border field
% 2008 science, https://www.science.org/doi/10.1126/science.1166466
% 1: firing rate higher than 0.3 times of the max
B_Score_v2=NaN(numCells,1);Cm=NaN(numCells,1);Dm=NaN(numCells,1);
 CC=[];tmp_Cm=[];tmp_Dm=[];tmpRM=[];tmpRM2=[];
for k=1:1:numCells
    tmpRM=GSrate_map{k}/max(GSrate_map{k},[],'all');
    tmpRM(tmpRM<0.3)=0;% to binary 
    tmpRM(tmpRM>=0.3)=1;% to binary
    
    % 2: Chose a peak with the max coverage of the given wall (Cm)
    bwl = bwlabel(tmpRM); %label binary object
    tmp_Cm=[];tmp_Cm2=[];
    for i=1:1: max(bwl,[],'all')
        if length(find(bwl==i))>20
            tmpRM2=zeros(50,50);
            tmpRM2(bwl==i)=1;

            % wall coverage
            tmp_Cm(i,1)=length(find(sum(tmpRM2(1:10,:),1) > 0));%north
            tmp_Cm(i,2)=length(find(sum(tmpRM2(40:50,:),1) > 0));%south
            tmp_Cm(i,3)=length(find(sum(tmpRM2(:,1:10),2) > 0)); %west
            tmp_Cm(i,4)=length(find(sum(tmpRM2(:,40:50),2) > 0));%east

            tmp_Cm2(i)=max(tmp_Cm(i,:))/50 ; %  coverage of the given wall
        else
            tmp_Cm(i,1)=NaN;%north
            tmp_Cm(i,2)=NaN;%south
            tmp_Cm(i,3)=NaN; %west
            tmp_Cm(i,4)=NaN;%east
            tmp_Cm2(i)=NaN ; %  coverage of the given wall
        end
    end
        
tmp_Dm=[];
    for i=1:1:max(bwl,[],'all')
        if length(find(bwl==i))>20
            tt=find(tmp_Cm(i,:)==max(tmp_Cm(i,:)));
            tt=tt(1);

            Weighted_RM=GSrate_map{k}/sum(GSrate_map{k}(bwl==i),'all');
            if tt==1 %north
                [X,Y] = meshgrid(1:1:50);
                weighted_Dmap=Y.*Weighted_RM/25;
            elseif tt==2 %south
                [X,Y] = meshgrid(50:-1:1);
                weighted_Dmap=Y.*Weighted_RM/25;
            elseif tt==3 %west
                [X,Y] = meshgrid(1:1:50);
                weighted_Dmap=X.*Weighted_RM/25;
            elseif tt==4 %east
                [X,Y] = meshgrid(50:-1:1);
                weighted_Dmap=X.*Weighted_RM/25;
            end

            tmp_Dm(i)=sum(weighted_Dmap(bwl==i)); 
            else
            tmp_Dm(i)=NaN;
        end

    end

    tmp_area=[];
    for i=1:1:max(bwl,[],'all')
         if length(find(bwl==i))>20
            tmp_area(i)=length(find(bwl==i));
         else
             tmp_area(i)=NaN;
         end
    end
    
    tmp_1=tmp_area.*( (tmp_Cm2-tmp_Dm)./ (tmp_Cm2 +tmp_Dm));
     B_Score_v2(k)= sum(tmp_1,'omitnan')/sum(tmp_area,'omitnan');

 
end

    Border_cells_v2=find(B_Score_v2>0.5);
    
save('ST_dF_grid_aut_data.mat','B_Score_v2','Border_cells_v2','-append');


%% visualize 
%{
for i=1:1:size(Border_cells_v2,1);
    k=Border_cells_v2(i);
    bs = round(B_Score_v2(k) , 2);
    
        subplot(2,2,1)
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',2);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 10);% Set 150ms delay
        daspect([vh vw vw])
        ylim([0 vh]);
        xlim([0 vw]);
         xticks([]); yticks([]);
        title({samplename;strcat("Cell#",  num2str(Original_Cell_ID(k)))},'FontSize',13);

        subplot(2,2,3) ;
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])
              mHz = round( max(GSrate_map{k},[],'all'), 1);
%         colorbar;
        title(strcat("Rate map, Hz=", num2str(mHz) ),'FontSize',13);

        subplot(2,2,[2 4])
        cla
        imshow(bg2,'Border','tight')
        hold on
        gscatter(CellPos(k,2)*s_dwSample,CellPos(k,1)*s_dwSample);
        title({ strcat("Border Score = ", num2str(bs) )},'FontSize',17);
        

        exportgraphics(gcf,strcat(CDIR,"\RateMap\", samplename,...
          " Cell#",num2str(Original_Cell_ID(k)),".jpg"),"Resolution",300)
        
    clf   
end 

CPmod=CellPos*s_dwSample;

CLR=B_Score_v2;
CLR(CLR<0)=0;CLR(find(isnan(CLR)))=0;
CLR=CLR/max(CLR);
clf
imshow(bg2,'Border','tight')
hold on
scatter(CPmod(:,2),CPmod(:,1),60,[1 1 1]);
for k=1:1:numCells
    scatter(CPmod(k,2),CPmod(k,1),60,[1 1 1-CLR(k)],'filled','MarkerFaceAlpha',CLR(k));
end

text(40,size(bg2,1)-40,strcat("Border cells, ", num2str(samplename)),'Color','white')

  exportgraphics(gcf,strcat(CDIR,"\Border Cell ", samplename,".jpg"))
close





for i=1:1:size(Border_cells,1);
    k=Border_cells(i);
    bs_1 = round(B_Score(k) , 2);
    bs_2 = round(B_Score_v2(k) , 2);
    

        
        subplot(2,2,1)
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',2);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 10);% Set 150ms delay
        daspect([vh vw vw])
        ylim([0 vh]);
        xlim([0 vw]);
         xticks([]); yticks([]);
        title({samplename;strcat("Cell#",  num2str(Original_Cell_ID(k)))},'FontSize',13);

        subplot(2,2,3) ;
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])
              mHz = round( max(GSrate_map{k},[],'all'), 1);
%         colorbar;
        title({strcat("Border Score (v1) = ", num2str(bs_1) );...
            strcat("Border Score (v2) = ", num2str(bs_2) );...
            strcat("Rate map, Hz=", num2str(mHz) )},'FontSize',13);


        exportgraphics(gcf,strcat(CDIR,"\RateMap_vsV1\", samplename,...
          " Cell#",num2str(Original_Cell_ID(k)),".jpg"))
    clf   
end 

%%
%}


for i=1:1:size(Border_cells_v2,1);
    k=Border_cells_v2(i);
    bs = round(B_Score_v2(k) , 2);
    
        subplot(2,1,1)
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',1);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 5);% Set 150ms delay
        daspect([vh vw vw])
        ylim([0 vh]);
        xlim([0 vw]);
         xticks([]); yticks([]);
        title({samplename;strcat("Cell#",  num2str(Original_Cell_ID(k)))},'FontSize',11);

        subplot(2,1,2) ;
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])
              mHz = round( max(GSrate_map{k},[],'all'), 1);
%         colorbar;
        title({strcat("BScore v2=", num2str(bs) );...
            strcat("Rate map, Hz=", num2str(mHz) )},'FontSize',11);

        

        exportgraphics(gcf,strcat(CDIR,"\Visualise\", samplename,...
          " Cell#",num2str(Original_Cell_ID(k)),".png"),"Resolution",600)
                exportgraphics(gcf,strcat(CDIR,"\Visualise\", samplename,...
          " Cell#",num2str(Original_Cell_ID(k)),".pdf"))
        
    clf   
end 


end




