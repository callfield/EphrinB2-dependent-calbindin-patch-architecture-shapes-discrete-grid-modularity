function Fn_meandF_Island_ex(DIR,ISLAND,samplename, s,t);

% 
% s=1
% t=1

% samplename=strcat(SampleName{3,s}," T",num2str(t))
% DIR=iwt_Dir{s,t};
% ISLAND=wt_Island{s,t};
dF= csvread(strcat(DIR, "\ST_PCI_noDup_dF.csv")); %
[NumFrame NumCell]=size(dF);
load(strcat(DIR,"\ST_dF_grid_aut_data.mat"), "caFr");
CDir=pwd;

ii=1;
for iii=unique(ISLAND).';
    Island{ii}=find(ISLAND==iii);
    ii=ii+1;
end
Island;
% normalize & movmean dF (max=1, min=0)
tmp_dF=dF-min(dF);tmp_dF=tmp_dF./max(tmp_dF);
tmp_dF=movmean(tmp_dF,10);

for i=1:3000:27001
    ST=i;
    ED=i+3000;
    dF_sum=tmp_dF(ST:ED,:).'+0.3;
ii=1;
    for iii=unique(ISLAND).';
        dF_i{ii}=tmp_dF(ST:ED,Island{ii}).'-0.3*iii;
        ii=ii+1;
    end

        clf
        hold on
        pSum=stdshade_sem(dF_sum,0.2,[0 0 0]);
ii=1;
    for iii=unique(ISLAND).';
        if iii==1;
            p1=stdshade_sem( dF_i{ii},0.2,[1 0 0]);
        elseif iii==2;
            p2=stdshade_sem( dF_i{ii},0.2,[0 1 1]);
        elseif iii==3;
            p3=stdshade_sem( dF_i{ii},0.2,[0 1 0]);
        elseif iii==4;
            p4=stdshade_sem( dF_i{ii},0.2,[1 0 1]);
        elseif iii==0;
            p0=stdshade_sem( dF_i{ii}-0.3*max(unique(ISLAND))-0.3,0.2,[0 0 0]);
        end
        ii=ii+1;
    end

    TICK= [0.3 -1*(unique(ISLAND).'*0.3+0.3)]+0.2;
    
    yticks([flip(TICK)]);
    % yticks([ mean(dF_p,"all" ) mean(dF_i4,"all" ) mean(dF_i3,"all" )   mean(dF_i2,"all" ) mean(dF_i1,"all" ) mean(dF_sum,"all" )  ]);
    if max(unique(ISLAND))==1
        yticklabels(["Other" "Island1" "All"]);
    elseif max(unique(ISLAND))==2
        yticklabels(["Other" "Island2" "Island1" "All"]);
    elseif max(unique(ISLAND))==3
        yticklabels(["Other" "Island3" "Island2" "Island1" "All"]);
    elseif max(unique(ISLAND))==4
        yticklabels(["Other" "Island4" "Island3" "Island2" "Island1" "All"]);
    end
        xticks([0:1000:3000]);xticklabels([((ST:1000:ED)-1)/caFr]);
        ax = gca; 
    ax.XAxis.FontSize = 12 ;
    ax.YAxis.FontSize = 14 ;
    xlim([0 3000])
    ylabel("mean dF(normalized)","FontSize",14)
    xlabel("Time(s)","FontSize",14)
    title({strcat(samplename)},'FontSize',15);
    exportgraphics(gcf,strcat(CDir, "\dFex_230609\dF_",samplename," f",num2str(ST),".jpg"));
    clf

end

end
