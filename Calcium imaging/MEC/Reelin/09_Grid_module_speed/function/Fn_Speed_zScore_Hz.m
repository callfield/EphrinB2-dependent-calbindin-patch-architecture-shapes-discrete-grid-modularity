function [ZSpeedHz, SpeedHz]=Fn_Speed_zScore_Hz(CDir)
close all;


load("ST_dF_grid_aut_data.mat")

x=split(pwd,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')))

caFr

%%
NumSes=length(sess);

Trk=[];All_dF=[];V=[];
for s= 1:1:NumSes; %
    strk= sTrk{s};
 
%     strk(:,4) = movmean(strk(:,4),caFr/2);
    tmp_V=strk(:,4);
    tmp_V(tmp_V==0)=NaN; % convert zero to nan before mean
    tmp_V = mean(tmp_V,round(caFr/2),'omitnan');

    Trk=[Trk;strk];
    V=[V;tmp_V];tmp_V=[];

   % dF and mean_dF
    dF{s} = (rawData3{s} - nanmean(rawData3{s}))./nanmean(rawData3{s}); % delta_d
%     Mean_dF{s} = movmean(dF{s} ,5); % moving average with upper 5 rows, row#1~4  
    All_dF = [All_dF;dF{s}];
%     AllMean_dF = [AllMean_dF;Mean_dF{s}];
    
end


All_dF(:,DupCell)=[];

%%
susp=find(Trk(:,2)==0&Trk(:,3)==0);%

Trk(susp,:)=[];Trk(:,1) =  1/caFr:1/caFr:size(Trk,1)/caFr;
All_dF(susp,:)=[];
V(susp,:)=[];
SMAX=round(max(V));
[numFrames numCells]=size(All_dF);
std_dF = nanstd(All_dF);
STD=2;

Vframe=cell(12,1);
for v=5:5:60
    Vframe{v/5}=find(V<=v & V>(v-5));
end
%%

% make reference z score using 0<v<=5
TestFrameNum=1000;
NumRepeat=10000; 
% prepare shuffle from 0<v<=5
RND=rand(TestFrameNum,NumRepeat);
% RND=fix( RND*length(Vframe{1}))+1;% only for z score using 0<v<=5
% RND=Vframe{1}(RND);% only for z score using 0<v<=5
RND=fix( RND*numFrames)+1; % z score using whole frame



% each cell
pk=cell(numCells,1);lk=cell(numCells,1);
p_from_varey=cell(numCells,1);w=cell(numCells,1);wxPk=cell(numCells,1);
FireFrame=cell(numCells,1);

SpeedHz=nan(12,numCells);
ZSpeedHz=nan(12,numCells);
for c=1:1:numCells 
    [pk{c},lk{c}, w{c}, p_from_varey{c},wxPk{c}] =  findpeaks_ho( All_dF(:,c), Trk(:,1),'MinPeakHeight', std_dF(:,c)*STD,'MinPeakProminence', std_dF(:,c));% 
    FireFrame{c} = round(lk{c}*caFr) ;

    tmp_Population=[];
    for rrr=1:NumRepeat

        tmp=intersect(FireFrame{c},RND(:,rrr));
        tmp_Population(rrr)=caFr*length(tmp)/TestFrameNum;
    end
%     histogram(tmp_Population)

%% 
    MEAN=nanmean(tmp_Population);
    SD=nanstd(tmp_Population);
    for v=5:5:60
        if length(Vframe{v/5})>=50;
            tmp=intersect(FireFrame{c},Vframe{v/5});
            SpeedHz(v/5,c)=caFr*length(tmp)/length(Vframe{v/5});
            ZSpeedHz(v/5,c)=(SpeedHz(v/5,c)-MEAN)/SD;
        end
    end

end


end



