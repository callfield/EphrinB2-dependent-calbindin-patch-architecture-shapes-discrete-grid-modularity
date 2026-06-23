clear all;close all

addpath(pwd);
addpath("function"));


load('../Data.mat');
CDir=pwd;
load(strcat(pwd,"\SpatialUniquness_Data.mat"));


CDir=pwd;



%%

zwt_mean=[];zwt_total=[];zwt_max=[];
n=1;
each_zwt_max=nan(2500,20);
for s=1:5
    for t=1:3
      
        zwt_mean=[zwt_mean; reshape(Mean_ErrDist_maxFRGrid{1,s,t},[],1)];
        zwt_total=[zwt_total; reshape(Total_ErrDist_maxFRGrid{1,s,t},[],1)];
        zwt_max=[zwt_max; reshape(Max_ErrDist_maxFRGrid{1,s,t},[],1)];
        if isempty(Mean_ErrDist_maxFRGrid{1,s,t})==0
            each_zwt_max(:,n)=reshape(Max_ErrDist_maxFRGrid{1,s,t},[],1);
            mean_zwt_max(:,n)=mean(reshape(Max_ErrDist_maxFRGrid{1,s,t},[],1),"omitmissing");
            name_wt(n)=strcat(SampleName{1,s}," t",num2str(t))
     
        n=n+1;
        end
    end
end

zeb_mean=[];zeb_total=[];zeb_max=[];
each_zeb_max=nan(2500,20);
n=1;
for s=1:7
    for t=1:3
      
        zeb_mean=[zeb_mean; reshape(Mean_ErrDist_maxFRGrid{2,s,t},[],1)];
        zeb_total=[zeb_total; reshape(Total_ErrDist_maxFRGrid{2,s,t},[],1)];
        zeb_max=[zeb_max; reshape(Max_ErrDist_maxFRGrid{2,s,t},[],1)];
        if isempty(Mean_ErrDist_maxFRGrid{2,s,t})==0
            each_zeb_max(:,n)=reshape(Max_ErrDist_maxFRGrid{2,s,t},[],1);
            mean_zeb_max(:,n)=mean(reshape(Max_ErrDist_maxFRGrid{2,s,t},[],1),"omitmissing");
             name_eb(n)=strcat(SampleName{2,s}," t",num2str(t))
     
        n=n+1;
        end

    end
end


 writematrix(zwt_max, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'WT_all',"AutoFitWidth",false)
 writematrix(zeb_max, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'EB2_all',"AutoFitWidth",false)
 writematrix(each_zwt_max, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'WT',"AutoFitWidth",false)
 writematrix(each_zeb_max, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'EB2',"AutoFitWidth",false)
 writematrix(mean_zeb_max, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'mean WT',"AutoFitWidth",false)
 writematrix(mean_zwt_max, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'mean EB2',"AutoFitWidth",false)
 writematrix(name_wt, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'WT name',"AutoFitWidth",false)
 writematrix(name_eb, strcat(CDir,'\GridMaxFR_MaxErroDist.xlsx'),'Sheet', 'EB2 name',"AutoFitWidth",false)



%%
clf
DATA1=log10(zwt_mean);
DATA2=log10(zeb_mean);
XLABEL="Mean Error Distance (log10(cm))";
SAMPLELABELS=["WT" "EB2-lacz"];
% YLABEL="Proportion";
TITLE="Mean Error distance (RateMap Corr Analysis)";
Fn_draw_cdf_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE)
exportgraphics(gcf,strcat(CDir, "\Mean Error Error Distance cdf WTvsEB2.jpg"),'Resolution',300)
exportgraphics(gcf,strcat(CDir, "\Mean Error Error Distance cdf WTvsEB2.pdf"))
clf
BIN=0.05;
Fn_draw_histogram_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE,BIN)
exportgraphics(gcf,strcat(CDir, "\Mean Error Error Distance hist WTvsEB2.jpg"),'Resolution',300)
exportgraphics(gcf,strcat(CDir, "\Mean Error Error Distance hist WTvsEB2.pdf"))
clf


%%



clf
DATA1=log10(zwt_total);
DATA2=log10(zeb_total);
XLABEL="Total Error Distance (log10(cm))";
SAMPLELABELS=["WT" "EB2-lacz"];
% YLABEL="Proportion";
TITLE="Total Error distance (RateMap Corr Analysis)";
Fn_draw_cdf_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE)
exportgraphics(gcf,strcat(CDir, "\Total Error Error Distance cdf WTvsEB2.jpg"),'Resolution',300)
exportgraphics(gcf,strcat(CDir, "\Total Error Error Distance cdf WTvsEB2.pdf"))
clf


BIN=0.05;
Fn_draw_histogram_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE,BIN)
exportgraphics(gcf,strcat(CDir, "\Total Error Error Distance hist WTvsEB2.jpg"),'Resolution',300)
exportgraphics(gcf,strcat(CDir, "\Total Error Error Distance hist WTvsEB2.pdf"))
clf


%%

clf
DATA1=log10(zwt_max);
DATA2=log10(zeb_max);
XLABEL="Max Error Distance (log10(cm))";
SAMPLELABELS=["WT" "EB2-lacz"];
% YLABEL="Proportion";
TITLE="Max Error distance (RateMap Corr Analysis)";
Fn_draw_cdf_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE)
exportgraphics(gcf,strcat(CDir, "\Max Error Error Distance cdf WTvsEB2.pdf"),'Resolution',300)

clf

% DATA1=zwt_max;
% DATA2=zeb_max;
% XLABEL="Max Error Distance(cm)";
BIN=0.05;
Fn_draw_histogram_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE,BIN)
exportgraphics(gcf,strcat(CDir, "\Max Error Error Distance hist WTvsEB2.jpg"),'Resolution',300)
clf

clf
DATA1=zwt_max;
DATA2=zeb_max;
XLABEL="Max Error Distance (cm)";
SAMPLELABELS=["WT" "EB2-lacz"];
% YLABEL="Proportion";
TITLE="Max Error distance (RateMap Corr Analysis)";
Fn_draw_cdf_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE)
exportgraphics(gcf,strcat(CDir, "\nolog Max Error Error Distance cdf WTvsEB2.jpg"),'Resolution',300)

clf

 % 
% DATA1=zwt_max;
% DATA2=zeb_max;
% XLABEL="Max Error Distance(cm)";
BIN=1;
Fn_draw_histogram_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE,BIN)
exportgraphics(gcf,strcat(CDir, "\nolog Max Error Error Distance hist WTvsEB2.jpg"),'Resolution',300)
clf
