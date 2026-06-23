function [all_relativeGridPhase, CSVlist_all, CSVlist]=AllPair_GridPhaseAnalysis(CDir)
close all;


x=split(pwd,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')));


%%

load('ST_dF_grid_aut_data.mat');

if length(Grid_Cells)>2



% fr = 10; % Frame rate for behavior video recording
vw = 100 ;% actual arena width(cm)
vh = 100 ;% actual arena height(cm)
% caFr = 10; % Frame rate for ca2+ recording
% vlim = 2; %defenition of moving
% nVoke1 or nVoke2
%Pfov_W = 1440;Pfov_H = 1080;Mfov_W = 900;Mfov_H = 600; %nVoke1
Pfov_W = 1280;Pfov_H = 800;Mfov_W = 950;Mfov_H = 600; %nVoke2
s_dwSample = 4 ; % rate of spatial downsample

bin = 2 ; %Analysis each (cm)
win = 5 ; % moving window size
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);

dF = csvread("ST_PCI_noDup_dF.csv"); %

[~, numCells] = size(dF);
CellPos = csvread("ST_PCI_noDup_CellPos.csv");
bg2=imresize(bg, s_dwSample, 'nearest');



%% all Pair Cross-correlation
% make pair of Grid Cell (cell with max grid score)

all_GridPair=[];
for k=1:1:length(Grid_Cells)
    for  l=1:1:length(Grid_Cells)
        tmp= [Grid_Cells(k),Grid_Cells(l)];
        all_GridPair=[all_GridPair;tmp];
    end
end



all_XCorr = cell(length(all_GridPair),1);
for k=1:1:length(all_GridPair) ;
    Ref=all_GridPair(k,1);
    Tar=all_GridPair(k,2);
    
% for k=1:1:5

  %  if length(m_lk{k}) > 100 ;
        for dx=-w_binNum:1:w_binNum
            for dy=-h_binNum:1:h_binNum
                aut1=zeros();aut2=zeros();aut3=zeros();aut4=zeros();aut5=zeros();aut6=zeros();aut7=zeros();
                n=0;
                   for x= 1:1:w_binNum
                       if 0<(x + dx) && (x + dx)*bin <=vw
                           for y=1:1:h_binNum
                              if 0<(y+dy) && (y+dy)*bin <= vh
                                  %Becareful!! (y,x) and y is start from upper left
                                  r1 = GSrate_map{Ref}(y, x);
                                  r2 = GSrate_map{Tar}(y+dy, x+dx);
                                aut1 = aut1 + r1*r2 ;
                                aut2 = aut2 + r1 ;
                                aut3 = aut3 + r2 ;
                                aut4 = aut4 + r1^2 ;
                                aut5 = aut5 + r1 ;
                                aut6 = aut6 + r2^2 ;
                                aut7 = aut7 + r2;    
                                n=n+1 ;
                              end    
                           end
                       end
                   end
                   if n>20
                       all_XCorr{k}(dy+h_binNum, dx+w_binNum) = (n*aut1-aut2*aut3)/(sqrt(n*aut4 - aut5^2) * sqrt(n*aut6-aut7^2)) ;                    
                   end

           end
  %      end

    end
end

save('ST_dF_grid_aut_data.mat','all_XCorr','all_GridPair','-append')






 all_relativeGridPhase=NaN(length(all_GridPair),9); 
 % 1,2:paired cell num, 3:Phase Distance, 4: Grid Scale Ratio(/Ref(1))
 % 5: physical distance
 % 6,7: ref/tar Grid scale, 
 % 8,9: ref/tar Grid score, 
 all_relativeGridPhase(:,1:2)=all_GridPair;
for k=1:1:length(all_GridPair) ;
    Ref=all_GridPair(k,1);
    Tar=all_GridPair(k,2);
    B = imgaussfilt(all_XCorr{k},3); % smoothed
    [cent, varargout]=FastPeakFind(B*255,0);
    [y x]=find(varargout==1);
    tmp_distance=sqrt((x*bin-w_binNum*bin).^2 + (y*bin-h_binNum*bin).^2);
    tmp_peak_id=knnsearch(tmp_distance,0,'K',1);
    all_relativeGridPhase(k,3)=tmp_distance(tmp_peak_id)/Grid_Scl_Ori(Ref,1);% Ration to reference grid scale
    all_relativeGridPhase(k,4)=Grid_Scl_Ori(Tar,1)/Grid_Scl_Ori(Ref,1);% Ration to reference grid scale
    
    
    all_relativeGridPhase(k,6)=Grid_Scl_Ori(Ref,1);% reference grid scale
    all_relativeGridPhase(k,7)=Grid_Scl_Ori(Tar,1);%target grid scale
    all_relativeGridPhase(k,8)=Grid_Score(Ref,1);% reference grid score
    all_relativeGridPhase(k,9)=Grid_Score(Tar,1);%target grid score
 
    
end

         if exist('CellPosBR', 'var') == 1
             cellDitance=sqrt((CellPosBR(all_GridPair(:,1),2)- CellPosBR(all_GridPair(:,2),2)).^2 + (CellPosBR(all_GridPair(:,1),1)-CellPosBR(all_GridPair(:,2),1)).^2 );
         else
             CellPosMod=zeros(numCells,2);
             CellPosMod(:,1)=CellPos(:,1)*s_dwSample*Mfov_W/Pfov_W;% DV position(um)
             CellPosMod(:,2)=CellPos(:,2)*s_dwSample*Mfov_H/Pfov_H;% ML position(um)
             cellDitance=sqrt((CellPosMod(all_GridPair(:,1),2)- CellPosMod(all_GridPair(:,2),2)).^2 + (CellPosMod(all_GridPair(:,1),1)-CellPosMod(all_GridPair(:,2),1)).^2 );
         end
         all_relativeGridPhase(:,5)=cellDitance;

save('ST_dF_grid_aut_data.mat','all_relativeGridPhase','-append')



         
                
         xx=0;i=1;list=[];
         CSVlist=NaN(length(cellDitance),1000/50);
         while xx <= max(cellDitance)
             tmp=cellDitance;
             tmp_id=find(tmp>xx & tmp< (xx+50) &all_relativeGridPhase(:,4)<1.1 & all_relativeGridPhase(:,4)>1/1.1 );
             if length(tmp_id)>1
                 list(i,1)=xx+25;
                 list(i,2)=mean(all_relativeGridPhase(tmp_id,3));
                 list(i,3)=std(all_relativeGridPhase(tmp_id,3));
                 % make Shuffle 
                 for RAND=1:1:100; %100 times Random sampling same number of pair  
                     rand_id=round((length(cellDitance)-1)*rand(length(tmp_id),1))+1;
                     rand_mean(RAND)=mean(all_relativeGridPhase(rand_id,3));% Shuffle
                 end
                 list(i,4)=mean(rand_mean);% Shuffle
                 list(i,5)=std(rand_mean);% Shuffle
             else
                 list(i,1)=xx+25;
                 list(i,2:5)=NaN;
                 
             end
             
             
             
             CSVlist(1:length(tmp_id),i)=all_relativeGridPhase(tmp_id,3);
            
             i=i+1;
             xx=xx+50;
         end
         
        
         
         errorbar(list(:,1),list(:,4),list(:,5));
         hold on
         errorbar(list(:,1),list(:,2),list(:,3),'k-d');
         legend("Shuffle","Real",'Location','northwest')
         
        xlabel("Physical Distance"); ylabel("Grid Phase Diff"); 
        title({samplename;strcat("AllPair Grid Phase Diff vs Cell Physical Distance (Scale ratio<1.1&>0.91)")});
       print(strcat(CDir,"\AllP_Phase_Diff\all Grid Phase Diff vs Physical Distance Similar Scale",samplename,".jpg"), '-djpeg', '-r0')
close


   
         xx=0;i=1;list=[];
         CSVlist=NaN(length(cellDitance),1000/50);
         while xx <= max(cellDitance)
             tmp=cellDitance;
             tmp_id=find(tmp>xx & tmp< (xx+50) );
             if length(tmp_id)>1
                 list(i,1)=xx+25;
                 list(i,2)=mean(all_relativeGridPhase(tmp_id,3));
                 list(i,3)=std(all_relativeGridPhase(tmp_id,3));
                 % make Shuffle 
                 for RAND=1:1:100; %100 times Random sampling same number of pair  
                     rand_id=round((length(cellDitance)-1)*rand(length(tmp_id),1))+1;
                     rand_mean(RAND)=mean(all_relativeGridPhase(rand_id,3));% Shuffle
                 end
                 list(i,4)=mean(rand_mean);% Shuffle
                 list(i,5)=std(rand_mean);% Shuffle
             else
                 list(i,1)=xx+25;
                 list(i,2:5)=NaN;
                 
             end
             
             
             
             CSVlist_all(1:length(tmp_id),i)=all_relativeGridPhase(tmp_id,3);
            
             i=i+1;
             xx=xx+50;
         end
         
       
         
         errorbar(list(:,1),list(:,4),list(:,5));
         hold on
         errorbar(list(:,1),list(:,2),list(:,3),'k-d');
         legend("Shuffle","Real",'Location','northwest')
         
        xlabel("Physical Distance"); ylabel("Grid Phase Diff"); 
        title({samplename;strcat("AllPair Grid Phase Diff vs Cell Physical Distance (All Grid)")});
       print(strcat(CDir,"\AllP_Phase_Diff\PW Grid Phase Diff vs Physical Distance All Grid",samplename,".jpg"), '-djpeg', '-r0')
close


else
    all_relativeGridPhase=[];
    CSVlist_all=[];
    CSVlist=[];
    
end

end
         
