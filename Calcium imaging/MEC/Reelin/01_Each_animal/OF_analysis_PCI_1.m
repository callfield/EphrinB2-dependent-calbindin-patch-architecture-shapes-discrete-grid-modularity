clear all
addpath("..\function\")
fr = 10; % Frame rate for behavior video recording
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
caFr = 10; % Frame rate for ca2+ recording
vlim = 2; %defenition of moving
Thresould=160;% To detect animal shape

% behavior video
vidDir = 'path\for\bahav\movie' ;
vidFileName{1} = 'bahavmovie.avi';


% behavior lamp video
vidDir2 = 'path\for\lamp\movie\' ;
LvidFileName{1} = 'lampmovie.avi';


vNum = length(LvidFileName)


%% tracking
trk = OF_videotrack_RevSusPreProcess(vidDir,vidFileName, vw, vh,Thresould);
save('RawTrack.mat','trk')
%load('PCI_grid_aut_data.mat')

%% Extracting Rec frame 
[Rec, RecStart] = OF_videotrack_lampPreProcess(vidDir2,LvidFileName, vw, vh);
save('RawTrack.mat','Rec','RecStart','vNum','-append') 



%% Rough visualise
mkdir Track
% separate to each session
trk2=[];
s=1;
for i=1:1:vNum
    [row, col]=size(RecStart{i});
   for k=1:1:row
       aaa= [RecStart{i}(k,1):RecStart{i}(k,2)];
       
       trk2= [trk2;trk{i}(aaa',:)];
       s=s+1;
   end
end
  susp=find(trk2(:,2)==0&trk2(:,3)==0&trk2(:,4)==0);%
  Trk=trk2  ;
  Trk(susp,:)=[];
  [numFrames, row] = size(Trk);
   Catime = numFrames/caFr ;
   Trk(:,1) = 1/caFr:1/caFr:Catime;


    plot(Trk(:,2),Trk(:,3),'k','LineWidth',2);
    ylim([0 vh]);
    xlim([0 vw]);
    daspect([vh vw vw])
    print(strcat("Track\all rec track.jpg"), '-djpeg', '-r0')
    close
    
    moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
    move_Trk = Trk(moveFrame,:);
    plot(move_Trk(:,2),move_Trk(:,3),'k','LineWidth',2);
    ylim([0 vh]);
    xlim([0 vw]);
    daspect([vh vw vw])
    print(strcat("Track\all moving rec track.jpg"), '-djpeg', '-r0')
    close
    
    h = plot(Trk(:,1), Trk(:,4), 'k');
    xlabel('time [sec]')
    ylabel('velocity [cm/sec]')
    print(strcat("Track\all rec velocity.jpg"), '-djpeg', '-r0')
    close