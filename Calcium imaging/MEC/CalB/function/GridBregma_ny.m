function [Slp_GS Slp_GW mod_CellPosDvM_SWO]= GridBregma_ny(SampleName,caFr,CV)

% zscore grid and 3 peaks
gscore_Z_patch;
G3_peak_patch;

%%

close all;
load('ST_dF_grid_aut_data.mat')


vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
vlim = 2; %defenition of moving
% STD=2; %Definition of peak
% %Pfov_W = 1440;Pfov_H = 1080;Mfov_W = 900;Mfov_H = 600; %nVoke1
Pfov_W = 1280;Pfov_H = 800;Mfov_W = 950;Mfov_H = 600; %nVoke2
s_dwSample = 4 ; % rate of spatial downsample
% 
bin = 2 ; %Analysis each (cm)
win = 5 ; % moving window size
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);

Trk = csvread("ST_PCI_Ca_behav_track.csv"); % 
dF = csvread("ST_PCI_noDup_dF.csv"); %
[numFrames numCells] = size(dF)
CellPos = csvread("ST_PCI_noDup_CellPos.csv"); % 1:DV, 2:ML


moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
move_dF = dF(moveFrame,:);
move_Trk = Trk(moveFrame,:);
std_dF = std(dF);




% Recal Grid Width
G_width=zeros(numCells,1);
for k=Grid_Cells
    [cent, varargout]=FastPeakFind(Autoc{k}*255*-1,0); ; %find local minima in autocorrecogram
    [row, col] = find(varargout==1);
    d=sqrt((row-h_binNum).^2 + (col-w_binNum).^2 );
    x=knnsearch(d,0,'K',1);
    G_width(k)=d(x)*2*bin; % distances from centerto 6 sorrounding peaks
end


mod_CellPosDvM_SWO=zeros(numCells,13);

mod_CellPosDvM_SWO(:,1)=CellPos(:,1)*1000/265+CV;% DV position(um)
mod_CellPosDvM_SWO(:,1)=1-1*mod_CellPosDvM_SWO(:,1)*cos(pi*12/180)-2050;% from Bregma, MEC tangential angle=12
h_ratio=Mfov_H/Pfov_H; 
w_ratio=Mfov_W/Pfov_W;
mod_CellPosDvM_SWO(:,2)=CellPos(:,2)*(h_ratio/ w_ratio)*(1000/265);% ML position(um)
mod_CellPosDvM_SWO(Grid_Cells,3:4)=Grid_Scl_Ori(Grid_Cells,:);% Grid Scale & Orientation
mod_CellPosDvM_SWO(Grid_Cells,5)= G_width(Grid_Cells,:);% Grid Width
for i=1:1:length(Grid_Cells)
    [B, I]=sort(G_3peaks{Grid_Cells(i)}(:,1)*180/pi);
    mod_CellPosDvM_SWO(Grid_Cells(i),6:8)=G_3peaks{Grid_Cells(i)}(I,1)*180/pi;% Three Orientation
    mod_CellPosDvM_SWO(Grid_Cells(i),9:11)=G_3peaks{Grid_Cells(i)}(I,2);% Grid scale
end
mod_CellPosDvM_SWO(:,12:13)=[Grid_Score Z_Grid_Score];
csvwrite("Grid_PosfromPAR_SWO.csv",mod_CellPosDvM_SWO)
save('ST_dF_grid_aut_data.mat','mod_CellPosDvM_SWO','-append');

gscatter(mod_CellPosDvM_SWO(Grid_Cells,3),mod_CellPosDvM_SWO(Grid_Cells,1))
h=lsline
Slp_GS=diff(h.YData)/diff(h.XData); % Slope of  least-squares line to Grid_Scale and DV position
close

gscatter(mod_CellPosDvM_SWO(Grid_Cells,5),mod_CellPosDvM_SWO(Grid_Cells,1))
h=lsline
Slp_GW=diff(h.YData)/diff(h.XData); % Slope of  least-squares line to Grid_Scale and DV position
close

CellPosBR(:,1:2)=mod_CellPosDvM_SWO(:,1:2);% DV position(um)


csvwrite("CellPosfromBregma_SWO.csv",CellPosBR)
save('ST_dF_grid_aut_data.mat','CellPosBR','-append');


end
