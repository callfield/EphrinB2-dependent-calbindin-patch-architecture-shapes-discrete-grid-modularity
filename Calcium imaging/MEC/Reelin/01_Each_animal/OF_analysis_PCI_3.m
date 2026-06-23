%%
clear all

load('ST_dF_grid_aut_data.mat')
addpath("..\function\")

fr = 10; % Frame rate for behavior video recording
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
caFr = 10; % Frame rate for ca2+ recording
STD=2; %Definition of peak
vlim = 2; %defenition of moving
% nVoke1 or nVoke2
%Pfov_W = 1440;Pfov_H = 1080;Mfov_W = 900;Mfov_H = 600; %nVoke1
Pfov_W = 1280;Pfov_H = 800;Mfov_W = 950;Mfov_H = 600; %nVoke2
s_dwSample = 4 ; % rate of spatial downsample

bin = 2 ; %Analysis each (cm)
win = 5 ; % moving window size
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);

Trk = csvread("ST_PCI_Ca_behav_track.csv"); % 
dF = csvread("ST_PCI_noDup_dF.csv"); %

[numFrames numCells] = size(dF)
CellPos = csvread("ST_PCI_noDup_CellPos.csv"); 


bg2=imresize(bg, s_dwSample, 'nearest');


SampleName='samplename';



%%  Detecting Ca2 + Fireing event during moving

mkdir ST_dF_Cell_activity ;


moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
move_dF = dF(moveFrame,:);
move_Trk = Trk(moveFrame,:);
std_dF = std(dF);

n = ceil(numFrames/9000);pk=cell(numCells,1);lk=cell(numCells,1);m_lk=cell(numCells,1);
p_from_varey=cell(numCells,1);w=cell(numCells,1);wxPk=cell(numCells,1);
for k=1:1:numCells 
%{
 for i =1:1:3 
        subplot(3,1,i)

        f_end = i*3000;
        f_str = f_end - 3000+1;        
        if i==n
            f_end = numFrames ;
        end
        x = Trk(f_str:f_end,1);
        y = dF(f_str:f_end,k); 
        [pks,locs] = findpeaks(y,x, 'MinPeakHeight', std_dF(:,k)*STD,'MinPeakProminence', std_dF(:,k));
              [pks,locs] = findpeaks(y,x, 'MinPeakHeight',std_Data_mean_d(:,k)*STD,'MinPeakProminence',std_Data_mean_d(:,k)*STD);
        plot(x,y,'k',locs,pks,'o') ;
        hline = refline([0 std_dF(:,k)*STD]);
        hline.Color = 'g';

end
   axes;
   title(strcat("Cell#",  num2str(Original_Cell_ID(k))));
   axis off;
   print(strcat("ST_dF_Cell_activity\Cell#",  num2str(Original_Cell_ID(k)), ".jpg"), '-djpeg', '-r0');
   close all 
  clf
%} 
    [pk{k},lk{k},w{k},p_from_varey{k},wxPk{k}] =  findpeaks_ho( dF(:,k), Trk(:,1),'MinPeakHeight', std_dF(:,k)*STD,'MinPeakProminence', std_dF(:,k));% save variables for prominence analysis

    fff = round(lk{k}*caFr) ;
    sp_v = Trk(fff,4) ; % velocity at firing 
    aaa = find(sp_v > 2) ; %velocity > 2cm/s
    m_lk{k} = round(lk{k}(aaa)*caFr) ; % frame of firing with velocity > 2cm/s
 %{
    clf
    plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',6); 
    ylim([0 vh]);
    xlim([0 vw]);
    hold on;
    plot(Trk( m_lk{k},2),Trk( m_lk{k},3),'r.', 'MarkerSize', 15)
    daspect([vh vw vw])
    title(strcat("Fireing map Cell#",  num2str(Original_Cell_ID(k)) ));
    print(strcat("ST_dF_Fireing_map\Fireing map Cell#",  num2str(Original_Cell_ID(k)), ".jpg"), '-djpeg', '-r0')
 %   close all
 %}
end
close all  

save('ST_dF_grid_aut_data.mat','lk', 'm_lk','-append');

%%  Fireing Rate Map during moving
x=move_Trk(:,2);y=move_Trk(:,3);
occup_map=zeros(w_binNum,h_binNum);
for i=1:1:w_binNum
    for l=1:1:h_binNum
        dd=sqrt(((x-i*bin).^2)+((y-l*bin).^2));
        occup_map(h_binNum-l+1,i)=length(find(dd<=win));% Count frame number within window 
     end
end
occup_map(occup_map==0) = 1 ;% convert 0 to 1
GSoccup_m = imgaussfilt(occup_map,2,'FilterSize', 5) ;

activ_m=cell(numCells,1);GSactiv_m=cell(numCells,1);GSrate_map=cell(numCells,1);

parfor k=1:1:numCells 
    x=Trk(m_lk{k},2);y=Trk(m_lk{k},3);
    activ_m{k}=zeros(w_binNum,h_binNum);;
    for i=1:1:w_binNum
        for l=1:1:h_binNum
            dd=sqrt(((x-i*bin).^2)+((y-l*bin).^2));
            activ_m{k}(h_binNum-l+1,i)=length(find(dd<=win));
        end
    end
    GSactiv_m{k} = imgaussfilt(activ_m{k},2,'FilterSize', 5);
    GSrate_map{k} = GSactiv_m{k}./(GSoccup_m/caFr);
end

save('ST_dF_grid_aut_data.mat','GSrate_map', 'GSoccup_m','-append');

disp("Rate map end")
%% Autocorrerogram

Autoc = cell(numCells,1);
for k=1:1:numCells ;
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
                                  r1 = GSrate_map{k}(y, x);
                                  r2 = GSrate_map{k}(y+dy, x+dx);
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
                       Autoc{k}(dy+h_binNum, dx+w_binNum) = (n*aut1-aut2*aut3)/(sqrt(n*aut4 - aut5^2) * sqrt(n*aut6-aut7^2)) ;                    
                   end

           end
  %      end

    end
end

save('ST_dF_grid_aut_data.mat','Autoc', '-append');
disp("Autocorrerogram end")
  
%% Grid score
%Transform Cartesian coordinates to polar or cylindrical
for k=1:1:numCells
    if isempty(Autoc{k})==0
        [w h]=size(Autoc{k});
        break
    end
end

[oy, ox] = ndgrid((w_binNum-1):-1:-(w_binNum-1),-(h_binNum-1):1:(h_binNum-1)); % make Cartesian coordinates, be careful
[theta,rho] = cart2pol(ox,oy);
tmp_rho = reshape(rho,[],1)*bin; % unit = bin
[unique_ids,~,idmatch_rho] = unique(tmp_rho); % aquire unique ID

Dis_Center_Autoc=cell(numCells,1);
in_r=cell(numCells,1);

    T_Autoc=cell(numCells,1);
    rot_T_Autoc=cell(numCells,6);
    Grid_Score=zeros(numCells,1);
    trim_Autoc=cell(numCells,1);
    rGS=zeros(numCells,round(vh));
    Cor_Autoc=cell(numCells,1);
for k=1:1:numCells
%    if length(m_lk{k}) > 80 ;  
    B = imgaussfilt(Autoc{k},3); % smoothed
    [cent, varargout]=FastPeakFind(B*255,0); %find peaks in autocorrecogram

    Angle=reshape(theta.*varargout,[],1);
    Dist=reshape(rho.*varargout,[],1);
    tmp=[Angle, Dist];
    G_all_peaks=tmp(tmp(:,1)~=0&tmp(:,2)~=0,:);
    G_6peaks=G_all_peaks(knnsearch(G_all_peaks(:,2),0,'K',6),:);
    in_r{k}=bin*min(G_6peaks(:,2))-10;
   

%Correct Elliptical distortions (Brandon et al., Science, 2011)
    G_3peaks=G_6peaks(G_6peaks(:,1)>0,:);
    r1=min(G_3peaks(:,2));
    R1=max(G_3peaks(:,2));
    min_Peak=G_3peaks(G_3peaks(:,2)==r1,:);
    max_Peak=G_3peaks(G_3peaks(:,2)==R1,:);
  % case1: minor axis=closet peak
    c = max_Peak(1,1)-min_Peak(1,1);%the shortest angle between the vector pointing to this field and the vector pointing to the farthest field

    Speculate_Major_axis=sqrt( (cos(c+pi/2)^2)/ (R1^-2-(sin(c+pi/2)/r1)^2));            
% compress along the major axis
    if Speculate_Major_axis/r1<2 &Speculate_Major_axis^2>0 % do not correct large distortion and % When peaks did not fit with ellipse
        phi_1=min_Peak(1,1)+pi/2; % 
        tilt_Autoc1=imrotate(Autoc{k},-phi_1*180/pi,'crop');% correct tilt
        T = [r1/Speculate_Major_axis    0  0; 0    1 0; 0    0  1];
        t_aff = affine2d(T);
        Cor_Autoc1 = imwarp(tilt_Autoc1,t_aff);  % compress along the major axis
    else
        Cor_Autoc1 = Autoc{k};
    end   
    

    
    [cor_y cor_x]=size(Cor_Autoc1); 
    [t_oy, t_ox] = ndgrid(-(cor_y-1)/2:1:(cor_y-1)/2,(cor_x-1)/2:-1:-(cor_x-1)/2); % make Cartesian coordinates, be careful
    [t_theta,t_rho] = cart2pol(t_ox,t_oy);
     B = imgaussfilt(Cor_Autoc1,3); % smoothed
    [cent, varargout]=FastPeakFind(B*255,0); %find peaks in autocorrecogram
     Angle=reshape(t_theta.*varargout,[],1);
    Dist=reshape(t_rho.*varargout,[],1);
    tmp=[Angle, Dist];
    G_all_peaks1=tmp(tmp(:,1)~=0&tmp(:,2)~=0,:);
    G_6peaks1=G_all_peaks1(knnsearch(G_all_peaks1(:,2),0,'K',6),:);
    in_r{k}=bin*min(G_6peaks1(:,2))-10;

% trim inner and outer mask 
    if isempty(in_r{k})==1 % when only one specific peak
        in_r{k}=30
    end
    in = t_rho*bin > round(in_r{k});
    out = t_rho*bin < round(in_r{k})+20;
    inout = in.*out ;
    trim_Autoc{k} = Cor_Autoc1.*inout;
    trim_Autoc{k} = trim_Autoc{k}./inout;
    r60= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},60,'crop'));
    r120= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},120,'crop'));
    r30= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},30,'crop'));
    r90= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},90,'crop'));
    r150= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},150,'crop'));

    tilt_GS1 = min([r60, r120]) - max([r30, r90, r150]);

% case2: major axis=farthest peak
    Speculate_Minor_axis=sqrt( (sin(c)^2)/ (1/r1^2-(cos(c)/R1)^2));


    % compress along the major axis
    if R1/Speculate_Minor_axis<2 &Speculate_Minor_axis^2>0 % do not correct large distortion
        phi_2=max_Peak(1,1); % 
        tilt_Autoc2=imrotate(Autoc{k},-phi_2*180/pi,'crop');% correct tilt
        T = [Speculate_Minor_axis/R1    0  0; 0    1 0; 0    0  1];
        t_aff = affine2d(T);
        Cor_Autoc2 = imwarp(tilt_Autoc2,t_aff);  % compress along the major axis
    else
        Cor_Autoc2 = Autoc{k};
    end   

    [cor_y cor_x]=size(Cor_Autoc2); 
    [t_oy, t_ox] = ndgrid(-(cor_y-1)/2:1:(cor_y-1)/2,(cor_x-1)/2:-1:-(cor_x-1)/2); % make Cartesian coordinates, be careful
    [t_theta,t_rho] = cart2pol(t_ox,t_oy);
     B = imgaussfilt(Cor_Autoc2,3); % smoothed
    [cent, varargout]=FastPeakFind(B*255,0); %find peaks in autocorrecogram
     Angle=reshape(t_theta.*varargout,[],1);
    Dist=reshape(t_rho.*varargout,[],1);
    tmp=[Angle, Dist];
    G_all_peaks2=tmp(tmp(:,1)~=0&tmp(:,2)~=0,:);
    G_6peaks2=G_all_peaks2(knnsearch(G_all_peaks2(:,2),0,'K',6),:);
    in_r{k}=bin*min(G_6peaks2(:,2))-10;
    
% trim inner and outer mask 
    if isempty(in_r{k})==1 % when only one specific peak
        in_r{k}=30
    end
    in = t_rho*bin > round(in_r{k});
    out = t_rho*bin < round(in_r{k})+20;
    inout = in.*out ;
    trim_Autoc{k} = Cor_Autoc2.*inout;
    trim_Autoc{k} = trim_Autoc{k}./inout;
    r60= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},60,'crop'));
    r120= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},120,'crop'));
    r30= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},30,'crop'));
    r90= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},90,'crop'));
    r150= nancorr2(trim_Autoc{k},  imrotate(trim_Autoc{k},150,'crop'));

    tilt_GS2 = min([r60, r120]) - max([r30, r90, r150]);
% chose bigger GS as Grid score
    rGS(k,round(in_r{k})+20) = max(tilt_GS1,tilt_GS2);
    if tilt_GS1>tilt_GS2
        Cor_Autoc{k}=Cor_Autoc1;
    else
        Cor_Autoc{k}=Cor_Autoc2;
    end

    if  isempty(max(rGS(k,rGS(k,:)~=0)))==0
        Grid_Score(k,1) = max(rGS(k,rGS(k,:)~=0));
    end

    Grid_Score(k,1) = max(rGS(k,rGS(k,:)~=0));

   %    end
end
save('ST_dF_grid_aut_data.mat','Grid_Score','Cor_Autoc','-append')
disp("Grid Score end")

hh=histogram(Grid_Score,'Normalization','probability','BinWidth',0.1);
x=hh.BinEdges(1:length(hh.BinEdges)-1);
plot(x,movmean(hh.Values,3),'k-','LineWidth',3)
xlim([-1 1.5])
xlabel('Gridness');
ylabel('Probability')
title(strcat(num2str(SampleName),", Grid score Distribution"));
print(strcat("Grid score Distribution.jpg"), '-djpeg', '-r0')
close

%% Caliculate significant Grid Score
RandNum=20 ;
[Gth, Rand_Grid_Score]= OF_grid_threshold_EllipticG3_noDup(RandNum, lk,vlim, Trk, GSoccup_m, caFr, numFrames, numCells,  vw, vh, bin, win,w,h);
save('ST_dF_grid_aut_data.mat','Gth','Rand_Grid_Score', '-append');  


RND=reshape(Rand_Grid_Score,[],1);
RND(RND==0)=NaN;
RND(find(isnan(RND)==1))=[];
Z_Grid_Score=(Grid_Score-mean(RND,'all'))/std(RND,0,'all');

save('ST_dF_grid_aut_data.mat','Z_Grid_Score','Grid_Score','-append')


%% visualize 
mkdir ST_dF_Autocorrerogram ;
% {
for k=1:1:numCells ;

        subplot(2,2,1)
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',2);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 10);% Set 150ms delay
        daspect([vh vw vw])
        ylim([0 vh]);
        xlim([0 vw]);
        title(strcat("Cell#",  num2str(Original_Cell_ID(k)) ));

        subplot(2,2,2) ;
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])
        colorbar;
        title(strcat("mov Ratio map Cell#",  num2str(Original_Cell_ID(k)) ));

        subplot(2,2,3)
        cla
        imshow(bg2,'Border','tight')
        hold on
        gscatter(CellPos(k,2)*s_dwSample,CellPos(k,1)*s_dwSample);
        title('Spatial Distribution of Grid cells in FoV');
        hold off


        subplot(2,2,4)
        cla
        imagesc(Autoc{k})
        xticks([]); yticks([]);
        daspect([1 1 1])
        colorbar;
        
        gs = round(Grid_Score(k) , 2);
        mHz = round( max(GSrate_map{k},[],'all'), 1);
       title(strcat("g= ", num2str(gs), ", ",num2str(mHz),"Hz", ", ","G_9_5= ", num2str(round(Gth,2))));
        print(strcat("ST_dF_Autocorrerogram\RateMapSummary Cell#",  num2str(Original_Cell_ID(k)), ".jpg"), '-djpeg', '-r0')

end    
       close all 
%}

%%

cellpos= [CellPos(:,2), CellPos(:,1)];
[hh ww rgb]=size(bg2);
SC=find(Grid_Score>Gth).';
SC=intersect(SC,1:1:numCells);
mkdir Grid_ST_dF
OUTDIR="Grid_ST_dF";

%% Prepair Cartesian coordinates to polar or cylindrical
for k=1:1:numCells
    if isempty(Autoc{k})==0
        [w h]=size(Autoc{k});
        break
    end
end
c_row = w_binNum-1;c_col = h_binNum-1;
[oy, ox] = ndgrid(c_col-(1:h)+1,(1:w)-c_row-1); % make Cartesian coordinates, be careful
[theta,rho] = cart2pol(ox,oy);% 


%% Define Grid cell using 6 nearst peaks
% Find peak in autocorrerogram
% FastPeak: https://www.mathworks.com/matlabcentral/fileexchange/37388-fast-2d-peak-finder
G_3peaks=cell(numCells,1);Grid_Scl_Ori=zeros(numCells,2);
NO_Grid=[]; 
for k=SC
    B = imgaussfilt(Autoc{k},3); % smoothed
    [cent, varargout]=FastPeakFind(B*255,0); ; %find peaks in autocorrecogram
    Angle=reshape(theta.*varargout,[],1);
    Dist=reshape(rho.*varargout,[],1);
    tmp=[Angle, Dist];
    G_all_peaks=tmp(find(tmp(:,1)~=0&tmp(:,2)~=0),:);
    G_6peaks=G_all_peaks(knnsearch(G_all_peaks(:,2),0,'K',6),:);
    
    tmp2=G_6peaks(knnsearch(G_6peaks(:,1),0,'K',4),1);% chose 4 peak
    Interaxis_Angle=abs( diff(sort(tmp2)));
    if min(Interaxis_Angle)<pi/6 || max(Interaxis_Angle)>pi/2 % if the minimal interaxis angle was smaller than 30�‹ or exceeded 90�‹
        NO_Grid=[NO_Grid;k];
        continue
    else
        tmp3=G_6peaks(knnsearch(G_6peaks(:,1),0,'K',3),2);% chose 3 peak
        diff_Gspace=tmp3./tmp3.';
        
        if isempty( find(diff_Gspace>2) ) ==0 % if far peaks are located too far
            NO_Grid=[NO_Grid;k];
            continue
        end
    end
    G_3peaks{k}=G_6peaks(knnsearch(G_6peaks(:,1),0,'K',3),:);% chose 3 peak

    Grid_Scl_Ori(k,1)=min(G_3peaks{k}(:,2))*bin; % Grid Space
    Grid_Scl_Ori(k,2)=mean(G_3peaks{k}(:,1))*180/pi; % Grid Orientation
    
    %{  
    [row, col] = find(varargout==1);
    d=sqrt((row-h_binNum).^2 + (col-w_binNum).^2 );
    x=knnsearch(d,0,'K',7);
    x(x==find(d==0))=[];
    G_6peaks(1:6,1)=row(x); % y coordinate of peaks to 6 sorrounding peaks
    G_6peaks(1:6,2)=col(x); % x coordinate of peaks to 6 sorrounding peaks 
    imagesc(Autoc{k}*255)
    hold on
    scatter(G_6peaks(:,2),G_6peaks(:,1),'k',"x")
    daspect([h_binNum,w_binNum,1])
    title(strcat("Surrounding 6 peaks, Cell",num2str(k)));
    print(strcat(num2str(OUTDIR),"\Cell ",num2str(k),".jpg"), '-djpeg', '-r0')
    close
    %}
end
    
Grid_Cells=setdiff(SC, NO_Grid);
x_scl=Grid_Scl_Ori(Grid_Scl_Ori(:,1)~=0,1);
y_ori=Grid_Scl_Ori(Grid_Scl_Ori(:,2)~=0,2);

save('ST_dF_grid_aut_data.mat','Grid_Cells','Grid_Scl_Ori','G_6peaks','-append')
 
 
for k=Grid_Cells
    cla
    B = imgaussfilt(Autoc{k},3); % smoothed
    [cent, varargout]=FastPeakFind(B*255,0); ; %find peaks in autocorrecogram
    Angle=reshape(theta.*varargout,[],1);
    Dist=reshape(rho.*varargout,[],1);
    tmp=[Angle, Dist];
    G_all_peaks=tmp(find(tmp(:,1)~=0&tmp(:,2)~=0),:);
    G_6peaks=G_all_peaks(knnsearch(G_all_peaks(:,2),0,'K',6),:);
    
    %{  
    [row, col] = find(varargout==1);
    d=sqrt((row-h_binNum).^2 + (col-w_binNum).^2 );
    x=knnsearch(d,0,'K',7);
    x(x==find(d==0))=[];
    G_6peaks(1:6,1)=row(x); % y coordinate of peaks to 6 sorrounding peaks
    G_6peaks(1:6,2)=col(x); % x coordinate of peaks to 6 sorrounding peaks 
    imagesc(Autoc{k}*255)
    hold on
    scatter(G_6peaks(:,2),G_6peaks(:,1),'k',"x")
    daspect([h_binNum,w_binNum,1])
    title(strcat("Surrounding 6 peaks, Cell",num2str(k)));
    print(strcat(num2str(OUTDIR),"\Cell ",num2str(k),".jpg"), '-djpeg', '-r0')
    
    %}
end
    close all


% subplot(4,1,[1:3])
% scatter(y_ori,x_scl,20,'filled')
% ylabel("Grid Scaling (cm)")
% xlabel("Orientation (deg)")
% title(strcat("Distribiution of Grid Scale(min)"));
%print(strcat(num2str(OUTDIR),"\Min Grid scale and orientation.jpg"), '-djpeg', '-r0')
%close
    
% Kernel smoothed density estimation
n=length(x_scl);

bandwid=2;
pdSix =fitdist(x_scl,'Kernel','BandWidth',bandwid);
x=linspace(min(x_scl)-20,max(x_scl)+20,100);
ySix = pdf(pdSix,x);
[peak locs]=findpeaks(ySix,x);

area(x,ySix,'FaceColor','k')
hold on
for i=1:1:length(locs)
    plot([locs(i);locs(i)],[min(ySix);max(ySix)*1.1] ,'c--','LineWidth',2)
    text(locs(i)-2,max(ySix)*1.15,num2str(round(locs(i),1)))
end
ylim([0 max(ySix)*1.2])
ylabel("Probability")
xlabel("Grid Scale (cm)")
title(strcat(num2str(SampleName),", KSD esrimate Grid Scale(Min)"));
print(strcat(num2str(OUTDIR),"\KSD Grid ScaleMean.jpg"), '-djpeg', '-r0')
close
 


% Kernel smoothed density estimation
n=length(y_ori);
bandwid=3.5;
pdSix =fitdist(y_ori,'Kernel','BandWidth',bandwid);
x=linspace(min(y_ori)-20,max(y_ori)+20,100);
ySix = pdf(pdSix,x);
[peak locs]=findpeaks(ySix,x);

area(x,ySix,'FaceColor','k')
hold on
for i=1:1:length(locs)
    plot([locs(i);locs(i)],[min(ySix);max(ySix)*1.1] ,'c--','LineWidth',2)
    text(locs(i)-2,max(ySix)*1.15,num2str(round(locs(i),1)))
end
ylim([0 max(ySix)*1.2])
ylabel("Probability")
xlabel("Grid Orientation (deg)")
title(strcat(num2str(SampleName),", KSD esrimate Grid Orientation"));
print(strcat(num2str(OUTDIR),"\KSD Grid Orientation.jpg"), '-djpeg', '-r0')
close
 


% summary
g=zeros(numCells,5);
for k=Grid_Cells
        g(k, :) = [k, cellpos(k,:), round(Grid_Scl_Ori(k,:))];
end
g=g(g(:,1)~=0,:);


subplot(6,12,1:6);
text(0,0.6,strcat('Number of detected cells =', num2str(length(1:1:numCells))))
text(0,0.3,strcat('Number of Grid cells =', num2str(length(Grid_Cells))))
axis off
subplot(6,12,7:12);
n=length(x_scl);
%a=10;b=105;
%bandwid=a/(n+b);
bandwid=2;
pdSix =fitdist(x_scl,'Kernel','BandWidth',bandwid);
x=linspace(min(x_scl)-20,max(x_scl)+20,100);
ySix = pdf(pdSix,x);
[peak locs]=findpeaks(ySix,x);

area(x,ySix,'FaceColor','k')
hold on
for i=1:1:length(locs)
    plot([locs(i);locs(i)],[min(ySix);max(ySix)*1.1] ,'c--','LineWidth',2)

end
xlim([min(g(:,4))-5 max(g(:,4))+5]);
ylabel("Probability")
xlabel("Grid Scale (cm)")
title(strcat(num2str(SampleName)));

 xticks([])

subplot(6,12,[13:18, 25:30, 37:42, 49:54, 61:66 ]);
imshow(bg2,'Border','tight')
hold on
gscatter(g(:,2)*s_dwSample,g(:,3)*s_dwSample,g(:,4));
text(g(:,2)*s_dwSample+1,g(:,3)*s_dwSample+2, num2str(g(:,4)),'Color','y','FontSize' ,8)
legend('off')

subplot(6,12,[19:24, 31:36, 43:48, 55:60, 67:72 ]);
gscatter(g(:,4),g(:,3)*s_dwSample);
ylim([0 hh])
xlim([min(g(:,4))-5 max(g(:,4))+5]);
axis ij
ylabel("Position in FoV (Dorsal to Ventral)")
xlabel("Minimum Grid Scale (cm)")
print(strcat(num2str(OUTDIR),"\FoV and grid scale (min) DV-axis.jpg"), '-djpeg', '-r0')
close


% summary(Grid Orientation)
subplot(6,12,1:6);
text(0,0.6,strcat('Number of detected cells =', num2str(length(1:1:numCells))))
text(0,0.3,strcat('Number of Grid cells =', num2str(length(Grid_Cells))))
axis off

subplot(6,12,7:12);
n=length(y_ori);
bandwid=3.5;
pdSix =fitdist(y_ori,'Kernel','BandWidth',bandwid);
x=linspace(min(y_ori)-20,max(y_ori)+20,100);
ySix = pdf(pdSix,x);
[peak locs]=findpeaks(ySix,x);

area(x,ySix,'FaceColor','k')
hold on
for i=1:1:length(locs)
    plot([locs(i);locs(i)],[min(ySix);max(ySix)*1.1] ,'c--','LineWidth',2)
end
ylim([0 max(ySix)*1.2])
xlim([min(g(:,5))-5 max(g(:,5))+5]);
ylabel("Probability")
xlabel("Grid Orientation (deg)")
title(strcat(num2str(SampleName)));
 xticks([])

subplot(6,12,[13:18, 25:30, 37:42, 49:54, 61:66 ]);
imshow(bg2,'Border','tight')
hold on
gscatter(g(:,2)*s_dwSample,g(:,3)*s_dwSample,g(:,5));
text(g(:,2)*s_dwSample+1,g(:,3)*s_dwSample+2, num2str(g(:,5)),'Color','y','FontSize' ,8)
legend('off')
subplot(6,12,[19:24, 31:36, 43:48, 55:60, 67:72 ]);
gscatter(g(:,5),g(:,3)*s_dwSample);
ylim([0 hh])
xlim([min(g(:,5))-5 max(g(:,5))+5]);
axis ij
ylabel("Position in FoV (Dorsal to Ventral)")
xlabel("Grid  Orientation (deg)")
print(strcat(num2str(OUTDIR),"\FoV and Grid Orientation.jpg"), '-djpeg', '-r0')
close


%% Grid Width
thr = 0.2; % Gwidth threshold

G_width = nan(numCells,1);
for c = 1:length(Grid_Cells)

    k = Grid_Cells(c);
    B = Autoc{k};

    % autocorrelogram center
    [nRow, nCol] = size(B);
    y0 = (nRow + 1) / 2;
    x0 = (nCol + 1) / 2;

    % distance from center for each pixel
    [X, Y] = meshgrid(1:nCol, 1:nRow);
    R = sqrt((X - x0).^2 + (Y - y0).^2);

    % nearest pixel where autocorrelation is below threshold
    maskBelow = isfinite(B) & B < thr & R > 0;

    if any(maskBelow(:))
        G_width(k) = min(R(maskBelow)) * 2 * bin;
    else
        G_width(k) = NaN;
    end
end




wid=zeros(numCells,5);
for k=Grid_Cells
    tmp=round(G_width(k));
    wid(k, :) = [k, cellpos(k,:), tmp, 1];
end
wid=wid(wid(:,1)~=0,:);



% Kernel smoothed density estimation
n=length(wid(:,4));
bandwid=3.5;
pdSix =fitdist(wid(:,4),'Kernel','BandWidth',bandwid);
x=linspace(min(wid(:,4))-20,max(wid(:,4))+20,100);
ySix = pdf(pdSix,x);
[peak locs]=findpeaks(ySix,x);

area(x,ySix,'FaceColor','k')
hold on
for i=1:1:length(locs)
    plot([locs(i);locs(i)],[min(ySix);max(ySix)*1.1] ,'c--','LineWidth',2)
    text(locs(i)-2,max(ySix)*1.15,num2str(round(locs(i),1)))
end
ylim([0 max(ySix)*1.2])
ylabel("Probability")
xlabel("Grid Width (cm)")
title(strcat(num2str(SampleName),", KSD esrimate Grid Width"));
print(strcat(num2str(OUTDIR),"\KSD Grid Width.jpg"), '-djpeg', '-r0')
close
 


% summary
subplot(6,12,1:6);
text(0,0.6,strcat('Number of detected cells =', num2str(length(1:1:numCells))))
text(0,0.3,strcat('Number of Grid cells =', num2str(length(Grid_Cells))))
axis off

subplot(6,12,7:12);
area(x,ySix,'FaceColor','k')
hold on
for i=1:1:length(locs)
    plot([locs(i);locs(i)],[min(ySix);max(ySix)*1.1] ,'c--','LineWidth',2)

end
ylim([0 max(ySix)*1.2])
xlim([min(wid(:,4))-5 max(wid(:,4))+5]);
ylabel("Probability")
xlabel("Grid Width (cm)")
title(strcat(num2str(SampleName)));
 xticks([])

subplot(6,12,[13:18, 25:30, 37:42, 49:54, 61:66 ]);
imshow(bg2,'Border','tight')
hold on
gscatter(wid(:,2)*s_dwSample,wid(:,3)*s_dwSample,wid(:,4));
text(wid(:,2)*s_dwSample+1,wid(:,3)*s_dwSample+3, num2str(wid(:,4)),'Color','y','FontSize' ,8)
legend('off')

subplot(6,12,[19:24, 31:36, 43:48, 55:60, 67:72 ]);
gscatter(wid(:,4),wid(:,3)*s_dwSample);
ylim([0 hh])
xlim([min(wid(:,4))-5 max(wid(:,4))+5]);
axis ij
ylabel("Position in FoV (Dorsal to Ventral)")
xlabel("Grid Width (cm)")
print(strcat(num2str(OUTDIR),"\FoV and grid width DV-axis.jpg"), '-djpeg', '-r0')
close




%visualise grid cell
for k=Grid_Cells ;
        subplot(2,2,1)
        cla
        plot(move_Trk(:,2),move_Trk(:,3),'Color',[0.8,0.8,0.8],'LineWidth',2);
        hold on;
        plot(Trk(m_lk{k},2),Trk(m_lk{k},3),'r.', 'MarkerSize', 10);% Set 150ms delay
        daspect([vh vw vw])
        ylim([0 vh]);
        xlim([0 vw]);
        title(strcat("Cell#",  num2str(Original_Cell_ID(k)) ));

        subplot(2,2,2) ;
        cla
        colormap(jet);
        imagesc(GSrate_map{k});
        xticks([]); yticks([]);
        daspect([vh vw vw])
        colorbar;
        title(strcat("mov Ratio map Cell#",  num2str(Original_Cell_ID(k)) ));

        subplot(2,2,3)
        cla
        imshow(bg2,'Border','tight')
        hold on
        gscatter(CellPos(k,2)*s_dwSample,CellPos(k,1)*s_dwSample);
        title('Spatial Distribution of Grid cells in FoV');
        hold off


        subplot(2,2,4)
        cla
        imagesc(Autoc{k})
        xticks([]); yticks([]);
        daspect([1 1 1])
        colorbar;
        
        gs = round(Grid_Score(k) , 2);
        mHz = round( max(GSrate_map{k},[],'all'), 1);

       title(strcat("g= ", num2str(gs), ", ",num2str(mHz),"Hz", ", ","G_9_5= ", num2str(round(Gth,2))));
        print(strcat(num2str(OUTDIR),"\RateMapSummary Cell#",  num2str(Original_Cell_ID(k)), ".jpg"), '-djpeg', '-r0')

end    
       close all 