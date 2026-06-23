function [Info_sum]=conditional_entropy_Zscore(CDir,caFr)
close all
% Info_sum
% 1: conditional entoropy (bit/sec)
% 2: conditional entoropy (bit/spike)
% 3: Normalized info
%% Load the data and variable.
x=split(pwd,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')))


Dir='';
load([Dir, 'ST_dF_grid_aut_data.mat'])

Trk = csvread([Dir, 'ST_PCI_Ca_behav_track.csv']); %
dF = csvread([Dir, 'ST_PCI_noDup_dF.csv']); %

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



%%
% nakashiba et al;
% https://www.science.org/doi/full/10.1126/science.1151120
% https://proceedings.neurips.cc/paper/1992/file/5dd9db5e033da9c6fb5ba83c7a7ebea9-Paper.pdf
% https://www.biorxiv.org/content/10.1101/2020.08.04.236174v3.full.pdf

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
I=NaN(numCells,1);I2=NaN(numCells,1);
Norm_I=NaN(numCells,1);
for k=1:1:numCells 
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
    
    prob_loc=GSoccup_m/sum(GSoccup_m,'all','omitnan');
    
    F=caFr*sum(GSactiv_m{k},'all','omitnan')/sum(GSoccup_m,'all','omitnan');

    H=prob_loc.*GSrate_map{k}.*log2(GSrate_map{k}/F);
    I(k)=sum(H,'all','omitnan'); % bit/sec
    I2(k)=I(k)/F;% bit/spike
    
    % make shuffle for 'normalized info' 
    % https://www.nature.com/articles/nature21692.pdf
    rand_I=NaN(100,1);
    for RAND=1:1:100
        rrr=ceil(rand(1,1)*numFrames);
        % ran_lk=ceil(rand(length(lk{k}),1)*numFrames);

        ran_lk=ceil(lk{k}*caFr)+rrr;
        ran_lk(find(ran_lk>numFrames))=ran_lk(find(ran_lk>numFrames))-numFrames;
        ran_lk=ceil(sort(ran_lk));

        aaa = Trk(ran_lk,4) ; % velocity at firing 
        bbb = find(aaa > vlim) ; %velocity > 2cm/s
        ran_m_lk = ran_lk(bbb) ; % frame of firing with velocity > 2cm/s
        
        ran_x=Trk(ran_m_lk,2);ran_y=Trk(ran_m_lk,3);
        ran_activ_m=zeros(w_binNum,h_binNum);
        for i=1:1:w_binNum
            for l=1:1:h_binNum
                dd=sqrt(((ran_x-i*bin).^2)+((ran_y-l*bin).^2));
                ran_activ_m(h_binNum-l+1,i)=length(find(dd<=win));
            end
        end
        ran_GSactiv_m = imgaussfilt(ran_activ_m,2,'FilterSize', 5);
        ran_GSrate_map = ran_GSactiv_m./(GSoccup_m/caFr);
        ran_F=caFr*sum(ran_GSactiv_m,'all','omitnan')/sum(GSoccup_m,'all','omitnan');

        ran_H=prob_loc.*ran_GSrate_map.*log2(ran_GSrate_map/ran_F);
        rand_I(RAND)=sum(ran_H,'all','omitnan'); % bit/sec
  
    end
    Norm_I(k)=(I(k)-mean(rand_I,'omitnan'))/std(rand_I,'omitnan');
  
end

    
Info_sum=[I, I2, Norm_I];
save([Dir, 'Spatial_Info.mat'],"Info_sum")


scatter(Z_Grid_Score,I2,'k','fill')
ylabel("Average spatial info (bit/spike)");xlabel("Grid Score(Zscore)");
title({samplename;"Average spatial info (bit/spike)"})
exportgraphics(gcf,strcat(CDir,"\conditional_entropy_zscore\",samplename," conditional info vs gridness.jpg"))
clf
scatter(Z_Grid_Score,Norm_I,'k','fill')
ylabel("Normalised Spatial Information");xlabel("Grid Score(Zscore)");
title({samplename;"Normalised Info=Conditional Entropy/mean(shuffled Conditional Entropy)"})
exportgraphics(gcf,strcat(CDir,"\conditional_entropy_zscore\",samplename," normalised info vs gridness.jpg"))
clf

 
 
end
