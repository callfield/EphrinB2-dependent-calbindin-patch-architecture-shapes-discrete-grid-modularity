function [Z_tmpGSrate_map]=Fn_zRatemap_forCorr_Mframe(DIR,mFRAME)

vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
bin = 2 ; %Analysis each (cm)
win = 5 ; % moving window size
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);


load(strcat(DIR, "\ST_dF_grid_aut_data.mat"),'m_lk','caFr')
numCells = size(m_lk,1);

vlim=2;
Trk = csvread(strcat(DIR, "\ST_PCI_Ca_behav_track.csv")); % 

moveFrame = find(Trk(:,4)>vlim); % velocity over limit (cm/s)
tmp_moveFrame=moveFrame(mFRAME);
move_Trk = Trk(tmp_moveFrame,:);


%%
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


activ_m=cell(numCells,1);GSactiv_m=cell(numCells,1);
tmpGSrate_map=nan(numCells,w_binNum,h_binNum);
for c=1:1:numCells 

    tmp_mlk=intersect(m_lk{c},tmp_moveFrame);

    x=Trk(tmp_mlk,2);y=Trk(tmp_mlk,3);
    activ_m{c}=zeros(w_binNum,h_binNum);;
    for i=1:1:w_binNum
        for l=1:1:h_binNum
            dd=sqrt(((x-i*bin).^2)+((y-l*bin).^2));
            activ_m{c}(h_binNum-l+1,i)=length(find(dd<=win));
        end
    end
    GSactiv_m{c} = imgaussfilt(activ_m{c},2,'FilterSize', 5);
    tmpGSrate_map_0 = GSactiv_m{c}./(GSoccup_m/caFr);
    tmpGSrate_map_0(GSoccup_m<10)= nan; % remove bin less than 10 frame
    tmpGSrate_map(c,1:w_binNum,1:h_binNum) = tmpGSrate_map_0;
end


Z_tmpGSrate_map=nan(numCells,50,50);
MEAN=mean(tmpGSrate_map,[2,3],"omitnan");
STD=std(tmpGSrate_map,[],[2 3],"omitnan");
Z_tmpGSrate_map=(tmpGSrate_map-MEAN)./STD;

end