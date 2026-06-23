function [SpeedCorr]=Speed_vs_meandFCorr2_NewP(CDir,caFr)
close all;


load("ST_dF_grid_aut_data.mat")

x=split(pwd,"\");
samplename=strcat(char(strrep(x(end-1),'_',' ')), " ",char(strrep(x(end),'_',' ')))

%%
NumSes=length(sess);

Trk=[];All_dF=[];ShiftV=[];
for s= 1:1:NumSes; %
    strk= sTrk{s};
 

    tmp_shiptF=[strk(:,4);repmat(0,10*5*2+1,1)];
    for t=1:1:10*5*2+1
        tmp=circshift(tmp_shiptF,round((t-10*5-1)*caFr/10));
       tmp_ShiftV(:,t)=tmp(1:end-(10*5*2+1));
    end
    Trk=[Trk;strk];
    ShiftV=[ShiftV;tmp_ShiftV];tmp_ShiftV=[];

   % dF and mean_dF
    dF{s} = (rawData3{s} - mean(rawData3{s}))./mean(rawData3{s}); % delta_d

    All_dF = [All_dF;dF{s}];

    
end


All_dF(:,DupCell)=[];

%%
susp=find(Trk(:,2)==0&Trk(:,3)==0);%

Trk(susp,:)=[];Trk(:,1) =  1/caFr:1/caFr:size(Trk,1)/caFr;
All_dF(susp,:)=[];
ShiftV(susp,:)=[];
SMAX=round(max(Trk(:,4)));
[numFrames numCells]=size(All_dF);
std_dF = nanstd(All_dF);
STD=2;


% each cell
pk=cell(numCells,1);lk=cell(numCells,1);
p_from_varey=cell(numCells,1);w=cell(numCells,1);wxPk=cell(numCells,1);
FireFrame=cell(numCells,1);FireTiming=zeros(numFrames,numCells);
vHZ=zeros(SMAX,numCells);
nonMoveFR=zeros(numCells,1);
for k=1:1:numCells 
    [pk{k},lk{k}, w{k}, p_from_varey{k},wxPk{k}] =  findpeaks_ho( All_dF(:,k), Trk(:,1),'MinPeakHeight', std_dF(:,k)*STD,'MinPeakProminence', std_dF(:,k));% 
    FireFrame{k} = round(lk{k}*caFr) ;


    for t=1:1:size(ShiftV,2)
            tlk= ShiftV(FireFrame{k},t);  
        for i=1:1:SMAX
            a=Trk(find(ShiftV(:,t)>i-1-4& ShiftV(:,t) <= i+4&ShiftV(:,t)>0),4);% moving window=10
            b=tlk(find(tlk>i-1-4 & tlk<=i+4 & tlk>0),1);% moving window=10

                      if length(a)<caFr*5
                    vHZ(i,k,t)=0;
                else
                    vHZ(i,k,t)=caFr*length(b)/length(a);
                    True_SMAX=i;
                end
        end
         [R(k,t) P(k,t)]=corr(vHZ(7:1:True_SMAX,k,t), (7:1:True_SMAX).','rows','complete');
    end
    
      a=Trk(find(ShiftV(:,51)<2 &ShiftV(:,51)>0),4);
      tlk= ShiftV(FireFrame{k},51);  
      b=tlk(find(tlk<2&tlk>0),1);% moving window=10  

     nonMoveFR(k)=caFr*length(b)/length(a);


end
%}

randNum=1000;

 rand_FireFrame=cell(numCells,1);
for k=1:1:numCells 
    tmp1=zeros(numFrames,1);
    tmp1(FireFrame{k})=1;
    rand_FireFrame{k}=zeros(length(FireFrame{k}),randNum);
    for RAND=1:1:randNum
        tmp2=circshift(tmp1,round(rand*numFrames));      
        rand_FireFrame{k}(:,RAND)=find(tmp2==1);
    end     
end

Pval=nan(numCells,1);
rand_R=nan(numCells,randNum);
for k=1:1:numCells 

    for RAND=1:1:randNum

        tlk= ShiftV(rand_FireFrame{k}(:,RAND));  
        for i=1:1:True_SMAX
        a=Trk(find(ShiftV(:,51)>i-1-4& ShiftV(:,51) <=i+4 & ShiftV(:,51)>0),4);% moving window=10
        b=tlk(find(tlk>i-1-4&tlk<=i+4&tlk>0),1);% moving window=10  
            if length(a)<caFr*5
                rand_vHZ(i,k)=0;
            else
                rand_vHZ(i,k)=caFr*length(b)/length(a);
                True_SMAX=i;
            end
        end
       [rand_R(k,RAND) rand_P]=corr(rand_vHZ(7:1:True_SMAX,k), (7:1:True_SMAX).','rows','complete');

    end

    if R(k,51)>mean(rand_R(k,:))
     Pval(k)=2*length(find( rand_R(k,:)>R(k,51) ))/randNum;
    elseif R(k,51)<mean(rand_R(k,:))
     Pval(k)=2*length(find( rand_R(k,:)<R(k,51) ))/randNum;
    end   
     

end
SpeedCorr{1}=R;
SpeedCorr{2}=P;
SpeedCorr{3}=vHZ;
SpeedCorr{4}=Pval;
SpeedCorr{5}=nonMoveFR;

clf
t=10*5+1;
rr=R(:,t);pp=P(:,t);
     k=find(rr==max(rr(Pval<=0.05))) ;k=k(1);
     
subplot(2,1,1)
    scatter(7:1:True_SMAX,vHZ(7:1:True_SMAX,k,t),25,'k','filled')
    xlabel("Velocity (cm/s)",'FontSize',15)   
    ylabel("Ca fireing rate (Hz)",'FontSize',15) 
%     [R(k) P(k)]=corr(y(3:1:SMAX), (3:1:SMAX).','rows','complete');
    FireNum=size(lk{k},1);
    HZ=round(FireNum*caFr/numFrames,3);
    title({strcat(samplename,", cell#",num2str(Original_Cell_ID(k)));...
        strcat("r = ",num2str(round(rr(k),2)),", p = ",num2str(Pval(k)),...
        ", Immobile FR=", num2str(round(nonMoveFR(k),2)),"Hz")},'FontSize',13)
subplot(2,1,2)
set(gca,'defaultAxesColorOrder',[[0 0.4470 0.7410]; [0 0 0]]);
yyaxis right
plot(1:100*caFr,Trk(1:100*caFr,4),'Color',[0.7 0.7 0.7],'LineWidth',2)
ylabel("Velocity",'FontSize',12)
hold on
yyaxis left
plot(1:100*caFr, All_dF(1:100*caFr,k),'LineWidth',3)
ylabel("dF/F",'FontSize',12)

xticks([100:100:100*caFr])
xticklabels([1:1:5])
xlabel("Time (s)",'FontSize',12)  

    exportgraphics(gcf,strcat(CDir,"\Speed_vs_Hz_NewP\SpeedVsFR ",samplename," max.jpg"),'Resolution',300)
    exportgraphics(gcf,strcat(CDir,"\Speed_vs_Hz_NewP\SpeedVsFR ",samplename," max.pdf"))
clf  

     k=find(rr==min(rr(Pval<=0.05))) ;
     subplot(2,1,1)
 scatter(7:1:True_SMAX,vHZ(7:1:True_SMAX,k,t),25,'k','filled')
    xlabel("Velocity (cm/s)",'FontSize',15)   
    ylabel("Ca fireing rate (Hz)",'FontSize',15) 
    FireNum=size(lk{k},1);
    HZ=round(FireNum*caFr/numFrames,3);
    title({strcat(samplename,", cell#",num2str(Original_Cell_ID(k)));...
        strcat("r = ",num2str(round(rr(k),2)),", p = ",num2str(Pval(k)),...
        ", Immobile FR=", num2str(round(nonMoveFR(k),2)),"Hz")},'FontSize',13)
    subplot(2,1,2)
set(gca,'defaultAxesColorOrder',[[0 0.4470 0.7410]; [0 0 0]]);
yyaxis right
plot(1:100*caFr,Trk(1:100*caFr,4),'Color',[0.7 0.7 0.7],'LineWidth',2)
ylabel("Velocity",'FontSize',12)
hold on
yyaxis left
plot(1:100*caFr, All_dF(1:100*caFr,k),'LineWidth',3)
ylabel("dF/F",'FontSize',12)

xticks([100:100:100*caFr])
xticklabels([1:1:5])
xlabel("Time (s)",'FontSize',12)    
    

    exportgraphics(gcf,strcat(CDir,"\Speed_vs_Hz_NewP\SpeedVsFR ",samplename," min.pdf"))
clf



end



