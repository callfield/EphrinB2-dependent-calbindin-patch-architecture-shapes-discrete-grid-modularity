function [filt_data]=GridPhase_vs_distance_sem_patch211028(CDir)
close all;

% Chose grid module for analysis

load("ST_dF_grid_aut_data.mat")

filt_data=[];cellDitance=[];
if length(Grid_Cells)>2
Trk = csvread("ST_PCI_Ca_behav_track.csv"); % 
CellPos = csvread("ST_PCI_noDup_CellPos.csv"); 
vlim = 2; %defenition of moving
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
Pfov_W = 1280;Pfov_H = 800;Mfov_W = 950;Mfov_H = 600; %nVoke2
moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
move_Trk = Trk(moveFrame,:);
s_dwSample = 4 ; % rate of spatial downsample
bg2=imresize(bg, s_dwSample, 'nearest');

% all_relativeGridPhase ( = AllP_GridPhase_WT, AllP_GridPhase_EB2)
 % 1,2:paired cell num, 3:Phase Distance, 4: Grid Scale Ratio(/Ref(1))
 % 5: physical distance
 % 6,7: ref/tar Grid scale, 
 % 8,9: ref/tar Grid score, 
 
 x=split(pwd,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')));

%%
% #1 similar Grid Scale,
Candidate=find(all_relativeGridPhase(:,4)<1.1 & ...
    all_relativeGridPhase(:,4)>0.91 );

filt_data=all_relativeGridPhase(Candidate,:);

% #2 remove grid cell its phase >0.7 (could not use for phase analysis)
A=filt_data(filt_data(:,3)>0.7, 1:2);
B = unique(A);

for i=1:1:length(B)  
    k=B(i);
    filt_data(filt_data(:,1)==k,:)=[];
    filt_data(filt_data(:,2)==k,:)=[];
end

% #3 ref has higher Zscore of Grid score 
   filt_data=filt_data(filt_data(:,8)> filt_data(:,9),:);


%%
if length(filt_data)>10
    xx=0;i=1;list=[];
    cellDitance=filt_data(:,5);
%          CSVlist=NaN(length(cellDitance),1000/50);
         while xx <= max(cellDitance)
             tmp=cellDitance;
             tmp_id=find(tmp>xx & tmp< (xx+30) );
             if length(tmp_id)>1
                 list(i,1)=xx+15;
                 list(i,2)=mean(filt_data(tmp_id,3));
                 list(i,3)=std(filt_data(tmp_id,3));
                 % make Shuffle 
                 for RAND=1:1:100; %100 times Random sampling same number of pair  
                     rand_id=round((length(cellDitance)-1)*rand(length(tmp_id),1))+1;
                     rand_mean(RAND)=mean(filt_data(rand_id,3));% Shuffle
                 end
                 list(i,4)=mean(rand_mean);% Shuffle
                 list(i,5)=std(rand_mean);% Shuffle
             else
                 list(i,1)=xx+15;
                 list(i,2:5)=NaN;
                 
             end

     
             i=i+1;
             xx=xx+30;
         end
         
      
         scatter(filt_data(:,5), filt_data(:,3),5,[0.8 0.8 0.8],'filled')
          hold on
         errorbar(list(:,1),list(:,4),list(:,5));
         errorbar(list(:,1),list(:,2),list(:,3),'k-d');
         legend("Real","Shuffle","Real (mean+sd)",'Location','northwest')

        xlabel("Physical Distance"); ylabel("Grid Phase Diff"); 
   title({samplename;strcat("Phase Diff vs Cell Physical Distance (Scale ratio<1.1&>0.91)")});
     exportgraphics(gcf,strcat(CDir,"\Phase_vs_distance\",samplename, "Grid Phase Diff vs Physical Distance Similar Scale.jpg"))

       close

end
    


tmp=filt_data(find(cellDitance>0 & cellDitance< 30),:);
tmp2=reshape(tmp(:,1:2).',[],1);

   for ii=1:1:length(tmp2)
       iii=(mod(ii,4)-1)*3+1;
       if mod(ii,4)==0
           iii=10;
       end


    k=tmp2(ii);

        subplot(13,12, [iii+12:iii+12+2 iii+12+12:iii+12+12+2 iii+12+24:iii+12+24+2  ] )
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',1);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 5);% Set 150ms delay
        daspect([vh vw vw])
        xticks([]); yticks([]);
        ylim([0 vh]); xlim([0 vw]);

        ZgScore = round(Z_Grid_Score(k) , 2);
        gScale = round(Grid_Scl_Ori(k) ,1);
        gPhase = round(tmp(ceil(ii/2),3),2);
        title({strcat("Cell#",  num2str(Original_Cell_ID(k)),", zg= ", num2str(ZgScore));
        strcat("gScale.= ",num2str(gScale));
        strcat("gPhase= ", num2str(gPhase))});

        subplot(13,12, [iii+48:iii+48+2 iii+48+12:iii+48+12+2 iii+48+24:iii+48+24+2  ] )
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])

        subplot(13,12, [iii+84:iii+84+2 iii+84+12:iii+84+12+2 iii+84+24:iii+84+24+2  ] )
        cla
        imagesc(Autoc{k})
        xticks([]); yticks([]);
        daspect([1 1 1])

        subplot(13,12, [iii+120:iii+120+2 iii+120+12:iii+120+12+2 iii+120+24:iii+120+24+2  ] )
        cla
        imshow(bg2,'Border','tight')
        hold on
        gscatter(CellPos(k,2)*s_dwSample,CellPos(k,1)*s_dwSample);
        
        if mod(ii,4)==0 || ii==length(tmp2)
              exportgraphics(gcf,strcat("Test2v2\Group", num2str(g)," high coef Cell#",  num2str(Original_Cell_ID(k)), ".jpg"))
            exportgraphics(gcf,strcat(CDir,"\Phase_vs_distance\",samplename, "closest cell pairs",  num2str(Original_Cell_ID(k)), ".jpg"))
             clf
            hold off
        end
   end
 
close

end
end




