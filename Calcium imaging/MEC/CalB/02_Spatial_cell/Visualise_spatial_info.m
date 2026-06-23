
function AllInfo=Visualise_spatial_info(CDir,DIR,SpatialInfo)
%% Load the data and variable.
x=split(DIR,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')))


Dir=DIR;
load([Dir, '\ST_dF_grid_aut_data.mat'])
load([Dir, '\Angle_FR_Score.mat'])
Trk = csvread([Dir, '\ST_PCI_Ca_behav_track.csv']); %
dF = csvread([Dir, '\ST_PCI_noDup_dF.csv']); %

vlim=2;
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
bin = 5 ; %Analysis each (cm)
win = 2.5 ; % moving window size
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);

% STD=2; % threshould to detect peak
[numFrames numCells]=size(dF);
move_Frames=find(Trk(:,4)>vlim);
move_Trk=Trk(move_Frames,:);
move_dF=dF(move_Frames,:);
[numMoveFrames numCells]=size(move_dF)




tmpID=find(SpatialInfo(:,3)>2);

for i=1:length(tmpID)
    k=tmpID(i);

        subplot(2,1, 1 )
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',1);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 5);% Set 150ms delay
        daspect([vh vw vw])
        xticks([]); yticks([]);
        ylim([0 vh]); xlim([0 vw]);

        ZgScore = round(Z_Grid_Score(k) , 2);
        Bscore = round(B_Score_v2(k) ,2);
        ninfo = round(SpatialInfo(k,3),2);
        title({strcat("Cell#",  num2str(Original_Cell_ID(k)),", zg= ", num2str(ZgScore));
        strcat("Bscore.= ",num2str(Bscore));
        strcat("Norm Info= ", num2str(ninfo))});

        subplot(2,1, 2 )
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])


            exportgraphics(gcf,strcat(CDir,"\allcell\",samplename, " Cell#",  num2str(Original_Cell_ID(k)), ".pdf"))
             exportgraphics(gcf,strcat(CDir,"\allcell\",samplename, " Cell#",  num2str(Original_Cell_ID(k)), ".jpg"),"Resolution",300)
             clf
            hold off

end

HD_S=nan(length(lk),1);
for c=1:length(lk)
   HD_S(c,1)=HD_Score_Bin6{c,1}(1,1);

end

AllInfo=[Z_Grid_Score, B_Score_v2, HD_S, SpatialInfo(:,3)];
end

