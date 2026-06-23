   close all; clear all

addpath('function');
addpath(pwd);


load('../08_GridModule/GridMod.mat'); 

load('..\Data.mat'); 

load('ZSpeedHz_each1.mat');

load('..\Data.mat');

CDir=pwd;


%%
mkdir ZScored_each1
gp1_ZSpeedHz_amean=cell(2,1);
gp1_ZSpeedHz_each=cell(2,1);
for s=1:5
    for t=1:3
        if isnan(gp1_GSratio(s,t)) ==0
            gcell= find(GROUP1{s,t}(:,4)~=0);

            for m=1:2
                tmp_ModGcell=gcell(gp1_Mod{s,t}{m});
    
                tmp=gp1_ZSpeedHz{s,t}(1:20,tmp_ModGcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean
                gp1_ZSpeedHz_amean{m}=[gp1_ZSpeedHz_amean{m}, nanmean(tmp,2)];
                gp1_ZSpeedHz_each{m}=[gp1_ZSpeedHz_each{m}, tmp];
            
            end
        end
    end
end

gp2_ZSpeedHz_amean=cell(2,1);
gp2_ZSpeedHz_each=cell(2,1);
for s=1:7
    for t=1:3
        if isnan(gp2_GSratio(s,t)) ==0

            gcell= find(GROUP2{s,t}(:,4)~=0);

            for m=1:2
                tmp_ModGcell=gcell(gp2_Mod{s,t}{m});
    
                tmp=gp2_ZSpeedHz{s,t}(1:20,tmp_ModGcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean
                gp2_ZSpeedHz_amean{m}=[gp2_ZSpeedHz_amean{m}, nanmean(tmp,2)];
                gp2_ZSpeedHz_each{m}=[gp2_ZSpeedHz_each{m}, tmp];
            
            end

        end
    end
end





writematrix(gp1_ZSpeedHz_each{1}, 'ZScored_each1_gp1_all.xlsx','Sheet', 'gp1_mod1',"AutoFitWidth",false)
writematrix(gp1_ZSpeedHz_each{2}, 'ZScored_each1_gp1_all.xlsx','Sheet', 'gp1_mod2',"AutoFitWidth",false)
writematrix(gp2_ZSpeedHz_each{1}, 'ZScored_each1_gp2_all.xlsx','Sheet', 'gp2_mod1',"AutoFitWidth",false)
writematrix(gp2_ZSpeedHz_each{2}, 'ZScored_each1_gp2_all.xlsx','Sheet', 'gp2_mod2',"AutoFitWidth",false)

%%

data1=gp1_ZSpeedHz_amean{1};




%%


clf
hold on
p1=stdshade_sem_nan2zero(gp1_ZSpeedHz_amean{1}.',0.2,[0 0 0]);
p2=stdshade_sem_nan2zero(gp1_ZSpeedHz_amean{2}.',0.2,[0 0 0.5]);

        ylabel("Z scored Frequency",'FontSize',15); 
        xlabel({"Speed (cm/s)"},'FontSize',20); 
        title(strcat("gp1, Grid cell"),'FontSize',20);
        xticks([0 5 10 15 20])
    xticklabels([0:5:20]);
      legend([p1 p2],{"Module 1", "Module 2"},'FontSize',15,'location','northwest');
        ax = gca;
        ax.XAxis.FontSize = 12;
        ax.YAxis.FontSize = 12;
         ylim([-1 2])
        % exportgraphics(gcf,strcat(CDir,"\GridMod_ZScored_DVspeed_gp1.jpg"))
        exportgraphics(gcf,strcat(CDir,"\GridMod_ZScored_DVspeed_gp1.pdf"))



clf
hold on
p1=stdshade_sem_nan2zero(gp2_ZSpeedHz_amean{1}.',0.2,[0 0 0]);
p2=stdshade_sem_nan2zero(gp2_ZSpeedHz_amean{2}.',0.2,[0.5 0 0]);

        ylabel("Z scored Frequency",'FontSize',15); 
        xlabel({"Speed (cm/s)"},'FontSize',20); 
        title(strcat("gp2-lz, Grid cell"),'FontSize',20);
xticks([0 5 10 15 20])
    xticklabels([0:5:20]);
      legend([p1 p2],{"Module 1", "Module 2"},'FontSize',15,'location','northwest')
        ax = gca;
        ax.XAxis.FontSize = 12;
        ax.YAxis.FontSize = 12;
        ylim([-1 2])
        % exportgraphics(gcf,strcat(CDir,"\GridMod_ZScored_DVspeed_gp2.jpg"))
        exportgraphics(gcf,strcat(CDir,"\GridMod_ZScored_DVspeed_gp2.pdf"))






         %%

writematrix(gp1_ZSpeedHz_amean{1}-gp1_ZSpeedHz_amean{2}, 'ZScored_Speed_DIF.xlsx','Sheet', 'gp1_mod1-2',"AutoFitWidth",false)
writematrix(gp2_ZSpeedHz_amean{1}-gp2_ZSpeedHz_amean{2}, 'ZScored_Speed_DIF.xlsx','Sheet', 'gp2_mod1-2',"AutoFitWidth",false)

