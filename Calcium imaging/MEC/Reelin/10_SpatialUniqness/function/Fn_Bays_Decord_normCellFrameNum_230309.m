function [ERedge_PN ERcent_PN ERall_PN]= Fn_Bays_Decord_normCellFrameNum_230309(DIR, num_AnaCellset, num_AnaFrameset)
close all


%% Load the data and variable.
t_step = 1; %temporal bin step for decording(s) 
num_AnaCellset=50;% number of cell in one analysis
num_AnaFrameset=10000;% number of frame for training in one analysis


X=split(DIR,"\");
samplename=strcat(char(strrep(X(end-1),'_',' ')), " ",char(strrep(X(end),'_',' ')))

load([DIR, '\ST_dF_grid_aut_data.mat'])

Trk = csvread([DIR, '\ST_PCI_Ca_behav_track.csv']); %
dF = csvread([DIR, '\ST_PCI_noDup_dF.csv']); %

vlim=2;
vw = 100 ;% acrutal arena width(cm)
vh = 100 ;% acrutal arena high(cm)
 bin = 10 ; %Analysis each (cm)
%  bin = 2 ; %Analysis each (cm)
w_binNum=ceil(vw/bin);
h_binNum=ceil(vh/bin);


STD=2; % threshould to detect peak
[numFrames numCells]=size(dF);
move_Frames=find(Trk(:,4)>vlim);
move_Trk=Trk(move_Frames,:);
move_dF=dF(move_Frames,:);
[numMoveFrames numCells]=size(move_dF)


if num_AnaCellset>numCells
    num_AnaCellset=numCells
end
if num_AnaFrameset>numMoveFrames
    num_AnaFrameset=numMoveFrames-1000
end

OutDir=strcat("Bayesean_Decode_AutoK_C",num2str(num_AnaCellset),"F",num2str(num_AnaFrameset));



%% Load the data.
%{
% %  spikes: 1xK cell array of spiking data, where K is the number of units:
% %         spikes{1}: vector of timestamps for each spike at unit 1.
% %         spikes{2}: vector of timestamps for each spike at unit 2.
% %         ...
% %   POS: NxM array of ground-truth stimulus features, where N is the number of
% %      stimulus dimensions and M is the number of sampling timestamps.
% %    t: 1xM array of timestamps for ground-truth stimulus sampling.
% %         sample_rate: sample rate of ground-truth data.
               >>  scale change from second to Hz (6 Hz has problem)

%}        
pk=cell(numCells,1);lk=cell(numCells,1);m_lk=cell(numCells,1);
p_from_varey=cell(numCells,1);w=cell(numCells,1);wxPk=cell(numCells,1);
for k=1:1:numCells 

    [pk{k},lk{k},w{k},p_from_varey{k},wxPk{k}] =  findpeaks_hisa_181225( dF(:,k), Trk(:,1),'MinPeakHeight',std(dF(:,k))*STD);% save variables for prominence analysis
    fff = round(lk{k}*caFr) ;
    zz=zeros(numFrames,1);
    zz(fff,1)=1;
    zz = zz(Trk(:,4)>vlim,1) ; % velocity > 2cm/s
    m_lk{k} = find(zz==1)/caFr ; % frame of firing with velocity > 2cm/s

end


allcell_spikes=m_lk;
POS=move_Trk(:,2:3).'; % "stimilous" means information intended to decode(like position)
t=(1:1:numMoveFrames)/caFr; % time stamp
sample_rate=caFr; % scale is Hz


%% Partition data for cross-validation(k-fold cross-validation)
% Fireing Rate Map during moving

% x=move_Trk(:,2);y=move_Trk(:,3);
% occup_map=zeros(w_binNum,h_binNum);occup_frame=cell(w_binNum,h_binNum);
% for i=1:1:w_binNum
%     for l=1:1:h_binNum
%         occup_frame{l,i}=find((i-1)*bin<x & x<i*bin &(l-1)*bin<y &y<l*bin);%  frame number within window 
%         occup_map(l,i)=length(find((i-1)*bin<x & x<i*bin &(l-1)*bin<y &y<l*bin));% Count frame number within window 
%      end
% end
% k:disjoint subsamples
% k=min(occup_map,[],'all'); % least number to include at least one binned location 
% if k<2
%     disp('Error: Need to inclease bin size')
%     k
% end
remFrame=numMoveFrames-num_AnaFrameset;
k=ceil(numMoveFrames/remFrame)
 
    


%%
for subgroup=1:1:k
    subgroup
    ErrD_PandN{subgroup}=[];
    
% Make subsample group (test group) of frame (divide by value of k)
    test_frame=(subgroup-1)*round(numMoveFrames/k)+1:1:subgroup*fix(numMoveFrames/k)-caFr*t_step;
    test_frame=test_frame(test_frame<=numMoveFrames);

%     test_frame_st=(test-1)*round(numMoveFrames/k)+1:caFr*t_step:test*round(numMoveFrames/k+caFr*t_step);
%     test_frame_st=test_frame_st(find(test_frame_st<=round(numMoveFrames-caFr*t_step)));

    tmp_train_use=repmat(1,numMoveFrames,1);
    tmp_train_use(test_frame,1)=0; 
    tmp2_train_use=find(tmp_train_use==1);
    trainlength=length(tmp2_train_use);
    
    for RAND=1:3 % training data starting from 3 different point

        train_frame=circshift(tmp2_train_use, round(trainlength/3),1);
        train_frame=sort(train_frame(1:num_AnaFrameset,1));
  %             train_frame=1:size(POS,2);      
        for cellset=1:1:ceil(numCells*3/num_AnaCellset)
 % select 50 cells randamly ten times
%             temp_cellset=randsample(1:1:numCells,num_AnaCellset);
             temp_cellset=Grid_Cells;
% temp_cellset=1:1:numCells;
            tmp_spikes=cell(num_AnaCellset,1);
            for s=1:1:length(temp_cellset)
                spikes{s}=allcell_spikes{temp_cellset(s)};
            end
        
%% Build the tuning curves.
            sigma = [3 3];
%             sigma =(10/bin)/4; % determine sigma based on filter size(11 cm). The default filter size of imgaussfilt is 2*ceil(2*sigma)+1.
            bin_size = [bin bin]; % Bin for Rate Map(X and Y axis, cm)
            f_base = 0.05;
            min_t_occ = 0.1;
%  [~,lambda] = build_tuning_curves_xValid(spikes, ... 
%     POS,t,sample_rate, ... 
%     train_frame, ... 
%     bin_size,sigma); 
             [~,lambda] = build_tuning_curves_xValid(spikes,POS,t,sample_rate,train_frame,bin_size,sigma);
%             [coords,alpha,beta] = build_tuning_curves_xValid(spikes,POS,t,sample_rate,train_frame,bin_size,sigma,f_base,min_t_occ);

%% Get information content curves.
            IC_curves = get_IC_curves(alpha,beta,f_base,min_t_occ);

%% Perform neural decoding.


            test_start = t(test_frame(1:2:end)); % Decode start from every 2 frame
            t_step=caFr*2;
            t_step=round(caFr/2);
            test_end = test_start + t_step;
%    test_start = t(test_frame(1)); % Decode start from every 2 frame
% t_step=caFr;
%             test_end = test_start + t_step;
            poiss_posterior = bayesian_decode_quantitate(spikes,test_start,test_end,lambda);

            [Z, X, Y] =size(poiss_posterior);
            poisPredPos=cell(Z,1);poisPredXY=zeros(Z,2);
            ActualXY=nan(Z,2);

binw=100/bin;
            StFrame=knnsearch((1:1:numMoveFrames).'/caFr.',test_start(1),'k',1);
            for z=1:1:Z
                data(1:binw,1:binw)=poiss_posterior(z,:,:);
                ActualXY(z,1:2)=POS(:,StFrame+2*z);% actual animal pos at test_start frame 
                tmpxy=ceil(ActualXY(z,1:2)/bin);
                PredProb(z)=data(tmpxy(1), tmpxy(2));
%          imshow(data,'InitialMagnification',5000);
                if isempty(find(data==max(data,[],'all')))==0
                    [poisPredXY(z,1), poisPredXY(z,2)]=find(data==max(data,[],'all'));
                else
                    poisPredXY(z,1)=NaN;
                end
            end


            poisPredXY=poisPredXY*bin-bin/2;% bin to actual pos
%             nbPredXY=nbPredXY*bin-bin/2;% bin to actual pos 

            ErrD_Poispred=sqrt((ActualXY(:,1)-poisPredXY(:,1)).^2+(ActualXY(:,2)-poisPredXY(:,2)).^2);
%             ErrD_Nbpred= sqrt((ActualXY(:,1)-nbPredXY(:,1)).^2+(ActualXY(:,2)-nbPredXY(:,2)).^2);
%             tmp=[ActualXY ErrD_Poispred ErrD_Nbpred];
%             ErrD_PandN{subgroup}=[ErrD_PandN{subgroup}; tmp];
        end
           

        
    end
end

%%

Err_occup_frame=cell(k,w_binNum,h_binNum);Err_AllDist=cell(k,w_binNum,h_binNum);
for subgroup=1:1:k
    X=ErrD_PandN{subgroup}(:,1);
    Y=ErrD_PandN{subgroup}(:,2);

    for z=1:1:w_binNum
         for x=1:1:h_binNum
            Err_occup_frame{subgroup,z,x}=find((z-1)*bin<X & X<z*bin &(x-1)*bin<Y &Y<x*bin);%  frame number within window     
            Err_AllDist{subgroup,z,x}=ErrD_PandN{subgroup}(Err_occup_frame{subgroup,z,x},3:4);
        end
    end
end

npErr=zeros(2,k,w_binNum,h_binNum);
npMinErr=cell(2);
npMinErr{1}=zeros(h_binNum,w_binNum);npMinErr{2}=zeros(h_binNum,w_binNum);
npMeanErr=cell(2);
npMeanErr{1}=zeros(h_binNum,w_binNum);npMeanErr{2}=zeros(h_binNum,w_binNum);

for z=1:1:w_binNum
    for x=1:1:h_binNum
        for subgroup=1:1:k
          total=length(Err_AllDist{subgroup,z,x}(:,1));
          p_err=length(find(Err_AllDist{subgroup,z,x}(:,1)>20));
          nb_err=length(find(Err_AllDist{subgroup,z,x}(:,2)>20));
          npErr(1,subgroup,z,x)=p_err/total ;
          npErr(2,subgroup,z,x)=nb_err/total ;
        end
        npMinErr{1}(h_binNum-x+1,z)=nanmin(npErr(1,:,z,x)) ;
        npMinErr{2}(h_binNum-x+1,z)=nanmin(npErr(2,:,z,x)) ;
        npMeanErr{1}(h_binNum-x+1,z)=nanmean(npErr(1,:,z,x)) ;
        npMeanErr{2}(h_binNum-x+1,z)=nanmean(npErr(2,:,z,x)) ;       
        
     end
end


% imshow(npMinErr{1},'InitialMagnification',10000)
% xlabel("x axis (cm)"); ylabel("y axis(cm)");
% c=colorbar
% title({strcat(SampleName);'Min error rate (Normalised cross-validated Bayesian decoding)' ;'Poisson encoding model'})
% c.Label.String = 'Error Rate';
% colormap('hot')
% print("Bayesean_Decode_AutoK\Normalised Min error rate map(Poisson).jpg", '-djpeg', '-r0')
% 
% clf
% imshow(npMinErr{2},'InitialMagnification',10000)
% xlabel("x axis (cm)"); ylabel("y axis(cm)");
% c=colorbar
% title({strcat(SampleName);'Min error rate (Normalised cross-validated Bayesian decoding)' ;'the negative binomial encoding model model'})
% c.Label.String = 'Error Rate';
% colormap('hot')
% print("Bayesean_Decode_AutoK\Normalised Min error rate map(NB model).jpg", '-djpeg', '-r0')

%% outpot
cd(CDir);
mkdir(OutDir);
save(strcat(OutDir,"\",samplename,' bays.mat'),'Err_AllDist','npMeanErr')


imshow(npMeanErr{1},'InitialMagnification',10000)
xlabel("x axis (cm)"); ylabel("y axis(cm)");
c=colorbar
title({samplename;strcat("Mean Error Rate (cell=",num2str(num_AnaCellset)," ,F=",num2str(num_AnaFrameset),')') ;'Bayesian decoding (Poisson model)'})
c.Label.String = 'Error Rate';
colormap('hot')
exportgraphics(gcf,strcat(OutDir,"\",samplename,"Mean ER map(Poisson).jpg"))

clf
imshow(npMeanErr{2},'InitialMagnification',10000)
xlabel("x axis (cm)"); ylabel("y axis(cm)");
c=colorbar
title({samplename;strcat("Mean Error Rate (cell=",num2str(num_AnaCellset)," ,F=",num2str(num_AnaFrameset),')') ;'Bayesian decoding (Negative Binomial model)'})
c.Label.String = 'Error Rate';
colormap('hot')
exportgraphics(gcf,strcat(OutDir,"\",samplename,"Mean ER map(NB model).jpg"))
close



% for Edge region csv
tmpEdge=zeros(size(npMeanErr{1}));
tmpEdge([1:20/bin w_binNum-20/bin+1:w_binNum],:)=1;
tmpEdge(:,[1:20/bin h_binNum-20/bin+1:h_binNum])=1;


% for center region csv
tmpCenter=zeros(size(npMeanErr{1}));
tmpCenter([(20/bin+1):(w_binNum-20/bin)],[(20/bin+1):(h_binNum-20/bin)])=1;

ERedge_PN=[npMeanErr{1}*tmpEdge, npMeanErr{2}*tmpEdge];
ERcent_PN=[npMeanErr{1}*tmpCenter, npMeanErr{2}*tmpCenter];
ERall_PN=[reshape(npMeanErr{1},[],1), reshape(npMeanErr{2},[],1)];

% csvwrite("edge_PoissonMeanErr.csv",npMeanErr{1}*tmpEdge); 
% csvwrite("edge_nbMeanErr.csv",npMeanErr{2}*tmpEdge); 
% csvwrite("center_PoissonMeanErr.csv",npMeanErr{1}*tmpCenter); 
% csvwrite("center_nbMeanErr.csv",npMeanErr{2}*tmpCenter); 
% csvwrite("all_PoissonMeanErr.csv",reshape(npMeanErr{1},[],1)); 
% csvwrite("all_nbMeanErr.csv",reshape(npMeanErr{2},[],1)); 


end
