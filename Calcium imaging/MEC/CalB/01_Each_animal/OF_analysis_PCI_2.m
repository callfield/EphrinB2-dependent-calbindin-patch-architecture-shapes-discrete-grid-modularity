
clear all
load('RawTrack.mat')

addpath("..\function\")


fr = 10; % Frame rate for behavior video recording
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
caFr = 10; % Frame rate for ca2+ recording
STD=3; %Definition of peak
vlim = 2; %defenition of moving
% nVoke1 or nVoke2
%Pfov_W = 1440;Pfov_H = 1080;Mfov_W = 900;Mfov_H = 600; %nVoke1
Pfov_W = 1280;Pfov_H = 800;Mfov_W = 950;Mfov_H = 600; %nVoke2
s_dwSample = 4 ; % rate of spatial downsample

bg = imread("Max_projectionof FOV.jpg");

filename = 'ST_dF_grid_aut_data.mat';

if exist(filename, 'file')
    save(filename, 'bg', '-append');
else
    save(filename, 'bg');
end
 
% Ca2+ csv data PCI
caDir = 'path/for/inscopix/output/csv' ;
caFileName = 'IDPS_output_raw_ca_trace.csv';
%Remove non-Cell signal
NOISE=["operationally removed cell"];
NOISE=NOISE+1;% IDPS start from cell#0


% Frame number in Each Ca Rec session 
sess(1:'num of ca rec session',1) = ["num frames of ca rec sessions"];
sess=sess+1;

NumSes = length(sess)
% cell position(PCI)
posFile = 'IDPS_output_-props.csv';
 


%% manually remove short recording



% separate to each session
trk2=[];sTrk= cell(1,NumSes); 
s=1;
for i=1:1:vNum
    [row, col]=size(RecStart{i});
   for k=1:1:row
       aaa= [RecStart{i}(k,1):RecStart{i}(k,2)];
       sTrk{s} = trk{i}(aaa',:);
       trk2=[trk2;trk{i}(aaa',:)];
       s=s+1;
   end
end


%% Read Ca2 signal
caNum = length(caFileName);
 Raw=[];
 Raw = csvread([caDir, char(caFileName)] ,2,1); % header, rowname �í�œ
 [numFrames, numCells] = size(Raw);
 Original_Cell_ID=1:1:numCells;
 Raw(:,NOISE)=[];%Remove Noise
 Original_Cell_ID(NOISE)=[]; 
 
 
 % Manually remove unused frame



% Separate to each context
rawData = cell(1,NumSes); dF= cell(1,NumSes); Mean_dF=  cell(1,NumSes);
pre_rawData=Raw;
for i =1:1:NumSes
    rawData{i} = pre_rawData(1:sess(i),:);
    pre_rawData(1:sess(i),:)=[];
    length(pre_rawData)
end

% Interporate dropped frame and shrink to [0 255] (remove and interporate NaN)
AllrawData = []; All_dF=[];AllMean_dF=[];rawData2 = cell(1,NumSes);
for s= 1:1:NumSes; %% chenge REC START
     [numFrames, numCells] = size(rawData{s});
     [vFrames, col] = size(sTrk{s});  
    aaa = 1:1:numFrames; vvv=linspace(1,numFrames,vFrames);
    rawData2{s}=zeros(vFrames, numCells);
    for i=1:1:numCells;
        ttt = 1:1:numFrames;
        nnn = find(isnan(rawData{s}(:,i)));
        ttt(nnn) = [];
        data = rawData{s}(:,i);
        data(nnn) = [];
     % Interporate dropped frame
        rawData{s}(:,i) = interp1(ttt, data, aaa, 'spline');
     % Change video length to fit with tracking record
        rawData2{s}(:,i) = interp1(aaa, rawData{s}(:,i), vvv, 'spline'); 
     % shrink to [0 255]   
        ddd = rawData2{s}(:,i) - min(rawData2{s}(:,i));
        rawData2{s}(:,i) = 255* ddd/max(ddd);

    end   

     % dF and mean_dF
    dF{s} = (rawData2{s} - mean(rawData2{s}))./ mean(rawData2{s}); % delta_d
    Mean_dF{s} = movmean(dF{s} ,5); % moving average with upper 5 rows, row#1~4 
    
  % Remove Suspecious Frame
    susp=find(sTrk{s}(:,2)==0&sTrk{s}(:,3)==0&sTrk{s}(:,4)==0);%
    rawData3{s}=rawData2{s};
    rawData2{s}(susp,:)=[];
    dF{s}(susp,:)=[];
    Mean_dF{s}(susp,:)=[];
    
     AllrawData = [AllrawData;rawData2{s}];
     All_dF=[All_dF;dF{s}];
     AllMean_dF = [AllMean_dF;Mean_dF{s}];
end
save('ST_dF_grid_aut_data.mat','rawData3','sTrk','RecStart','sess','-append')


Mean_dF=AllMean_dF;dF=All_dF;
csvwrite("ST_PCI_dF.csv", dF);



%% Cell location
CellPos = zeros(numCells,2);
Raw = csvread([caDir, char(posFile)] ,1,2); % header, rowname �í�œ
CellPos =  Raw(:,4:5);
CellPos(NOISE,:)=[];%Remove Noise

%% Make Track file
Trk=[];
for s= 1:1:NumSes; %% chenge RECSTART
    strk= sTrk{s};
    susp=find(sTrk{s}(:,2)==0&sTrk{s}(:,3)==0&sTrk{s}(:,4)==0);%
    strk(susp,:)=[];
    Trk=[Trk;strk];
end
[numFrames, numCells] = size(Mean_dF);
Catime = numFrames/caFr ;
Trk(:,1) = 1/caFr:1/caFr:Catime;
csvwrite("ST_PCI_Ca_behav_track.csv", Trk);


% visualise
%{
    mkdir Track
    plot(Trk(:,2),Trk(:,3),'k','LineWidth',2);
    ylim([0 vh]);
    xlim([0 vw]);
    daspect([vh vw vw])
    print(strcat("Track\all track.jpg"), '-djpeg', '-r0')
    close
    
    moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
    move_Trk = Trk(moveFrame,:);
    plot(move_Trk(:,2),move_Trk(:,3),'k','LineWidth',2);
    ylim([0 vh]);
    xlim([0 vw]);
    daspect([vh vw vw])
    print(strcat("Track\all moving track.jpg"), '-djpeg', '-r0')
    close
%}
%% Remove same cell signal (use Mean_dF)

% Thresould 1: distance between cells (within 6pics)

D=pdist(CellPos);
[row,col] =find(squareform(D)<6 & squareform(D)>0); % omit same cell

% Thresould 2: high correlation (r>0.5) of dF

if isempty(row)==0
   
    r=zeros(length(row),1);
   for i=1:1:length(row)
    r(i)=corr(dF(:,row(i)),dF(:,col(i)));
   end
   xx = find(r>0.5);
end


% Thresould 3: keep a cell with higher max signal
   c=[];
    for i=1:1:length(xx)
        a = row(xx(i));
        b = col(xx(i));

        MaxA = max(dF(:,a));
        MaxB = max(dF(:,b));
        if MaxA > MaxB
            c = [c;b];
        else
            c = [c;a] ;
        end
    end
    DupCell=unique(c);
    % remove duplicated cell
     Cells=1:1:numCells;
   Cells(DupCell)=[];
   
csvwrite("ST_PCI_noDup_dF.csv", dF(:,Cells));
csvwrite("ST_PCI_noDup_CellPos.csv",CellPos(Cells,:)); 



% visualise
bg2=imresize(bg, 2*s_dwSample, 'nearest');
imshow(bg2,'Border','tight')
hold on
%scatter(CellPos(:,1)*s_dwSample,CellPos(:,2)*s_dwSample,15,'c','filled');
text(CellPos(:,2)*2*s_dwSample,CellPos(:,1)*2*s_dwSample,num2str(Original_Cell_ID.'),'Color','cyan')
hold on
scatter(CellPos(row(xx),2)*2*s_dwSample,CellPos(row(xx),1)*2*s_dwSample,35,'y');
 hold on
scatter(CellPos(DupCell,2)*2*s_dwSample,CellPos(DupCell,1)*2*s_dwSample,35,'r');

print(strcat("Removed Cells.jpg"), '-djpeg', '-r0')
close


Original_Cell_ID(DupCell)=[];
csvwrite("Original_Cell_ID_noDup.csv",Original_Cell_ID); 
save('ST_dF_grid_aut_data.mat', 'Cells','Original_Cell_ID','DupCell','-append')

