close all; clear all

load('../Data.mat');

addpath(pwd)


CDir=pwd;


mkdir conditional_entropy_zscore


Spatial_Info_Group1=cell(7,3);
Spatial_Info_Group2=cell(7,3);

% Spatial_Info_Group1 & Spatial_Info_Group2
% 1: confitional entoropy (bit/sec)
% 2: confitional entoropy (bit/spike)
% 3: Normalised info (bit/sec)
%%
 caFr=10;
for s=1:7
    for t=1:3
        cd(group1_Dir{s,t});
        [Spatial_Info_Group1{s,t}]=conditional_entropy_Zscore(CDir,caFr);

    end
end


for s=1:7
    for t=1:3

        cd(group2_Dir{s,t});
        [Spatial_Info_Group2{s,t}]=conditional_entropy_Zscore(CDir,caFr);
                                                                

    end
end

save(strcat(CDir,"/Spatial_Info.mat"),"Spatial_Info_Group2","Spatial_Info_Group1")


%% 
gru1=[];
for s=1:7
    for t=1:3
        tmp=Spatial_Info_Group1{s,t};
        gru1=[gru1;tmp];
    end
end

gru2=[];
for s=1:7
    for t=1:3

        tmp=Spatial_Info_Group2{s,t};
        gru2=[gru2;tmp];

    end
end


writematrix(gru1(:,3),"NormSpatial_Info.xlsx","Sheet","GROUP1")
writematrix(gru2(:,3),"NormSpatial_Info.xlsx","Sheet","GROUP2")


%% visualise spatial info (bits/spike)
    tmp_gru1=gru1(:,2);
    tmp_gru2=gru2(:,2);
    
h1=cdfplot(tmp_gru1)
hold on
h2=cdfplot(tmp_gru2)
h1.LineWidth=2;h2.LineWidth=2
h1.Color="black"

xlabel("Average spatial info (bits/spike)")

ylabel("Proportion");
title("")
 yticks([0 0.5 1]);
box off

hold on
[f1,x1] = ecdf(tmp_gru1);
[f2,x2] = ecdf(tmp_gru2);
Idx = knnsearch(x2,x1);
ksd=abs(f2(Idx)-f1);
ksd_id=find(ksd==max(ksd));
plot([x1(ksd_id) x2(Idx(ksd_id))],[f1(ksd_id) f2(Idx(ksd_id))],'k--','LineWidth',1)

grid off
legend("GROUP1","GROUP2","location","northwest")

[h,p,ks2stat] = kstest2(tmp_gru1, tmp_gru2)
title({"Average spatial info (bit/spike)";strcat("p = ",num2str(p), ", KS distance = ", num2str(ks2stat) )})
    exportgraphics(gcf,strcat(CDir,"\conditional_entropy_zscore\Average spatial info per spike cdf GROUP1 vs GROUP2.pdf"))

close
  

%% Normalised info (bit/sec)
    tmp_gru1=gru1(:,3);
    tmp_gru2=gru2(:,3);
    
h1=cdfplot(tmp_gru1)
hold on
h2=cdfplot(tmp_gru2)
h1.LineWidth=2;h2.LineWidth=2
h1.Color="black"

xlabel("Normalised info (bit/sec)")

ylabel("Proportion");
title("")
 yticks([0 0.5 1]);
box off

hold on
[f1,x1] = ecdf(tmp_gru1);
[f2,x2] = ecdf(tmp_gru2);
Idx = knnsearch(x2,x1);
ksd=abs(f2(Idx)-f1);
ksd_id=find(ksd==max(ksd));
plot([x1(ksd_id) x2(Idx(ksd_id))],[f1(ksd_id) f2(Idx(ksd_id))],'k--','LineWidth',1)

grid off
legend("GROUP1","GROUP2","location","northwest")

[h,p,ks2stat] = kstest2(tmp_gru1, tmp_gru2)
title({"Normalized spatial info.";strcat("p = ",num2str(p), ", KS distance = ", num2str(ks2stat) )})
    exportgraphics(gcf,strcat(CDir,"\conditional_entropy_zscore\Normalized spatial info cdf GROUP1 vs GROUP2.pdf"))

close
  



