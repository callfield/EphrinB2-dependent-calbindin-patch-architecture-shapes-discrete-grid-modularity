    clear all;close all


addpath(pwd);

addpath('function');
load('../Data.mat');

CDir=pwd;


%% output: gp1_ZSpeedHz
% Zscored based on Hz on v=0-5 cm/s (Fn_Speed_zScore_Hz)
% gp1_SpeedHz
% Hz on each 5cm/s (Fn_Speed_zScore_Hz)
%%
mkdir ZScored_all_trial


gp1_ZSpeedHz=cell(5,3);gp1_SpeedHz=cell(5,3);
for s=1:5
    for t=1:3
        samplename=strcat(SampleName{1,s}, " T",num2str(t));
         cd(group1_Dir{s,t})   
         [gp1_ZSpeedHz{s,t}, gp1_SpeedHz{s,t}]=Fn_Speed_zScore_Hz(CDir);


        for dv=1:1:3
            STR=200*dv+1850;
            END=200*(dv+1)+1850;
            tmpcell= find(GROUP1{s,t}(:,1)<-1*STR & GROUP1{s,t}(:,1)>-1*END);
            if isempty(tmpcell)==0

                tmp=gp1_ZSpeedHz{s,t}(1:4,tmpcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean

                stdshade_sem_nan2zero(tmp.',0.2,...
                    [0+(dv-1)*0.5 0+(dv-1)*0.05 0+(dv-1)*0.05])
                hold on
            end

        end
        
           

        ylabel("Zscored Frequency",'FontSize',15); 
        xlabel({"Speed (cm/s)"},'FontSize',20); 
        title(strcat(samplename,", All cell"),'FontSize',20);
        xticklabels(["0-5" "5-10" "10-15" "15-20" ])
        ax = gca;
        ax.XAxis.FontSize = 12
        ax.YAxis.FontSize = 12
        exportgraphics(gcf,strcat(CDir,"\ZScored_all_trial\ZScored_DVspeed_",samplename,".jpg"))
clf
    end
end


save("ZSpeedHz.mat","gp1_ZSpeedHz","gp1_SpeedHz",'-append')

gp2_ZSpeedHz=cell(7,3);gp2_SpeedHz=cell(7,3);
for s=1:7
    for t=1:3
        samplename=strcat(SampleName{2,s}, " T",num2str(t));

         cd(group2_Dir{s,t})   
        [gp2_ZSpeedHz{s,t}, gp2_SpeedHz{s,t}]=Fn_Speed_zScore_Hz(CDir);

        for dv=1:1:3
            STR=200*dv+1850;
            END=200*(dv+1)+1850;
            tmpcell= find(GROUP2{s,t}(:,1)<-1*STR & GROUP2{s,t}(:,1)>-1*END);
            if isempty(tmpcell)==0

                tmp=gp2_ZSpeedHz{s,t}(1:4,tmpcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean

                stdshade_sem_nan2zero(tmp.',0.2,...
                    [0+(dv-1)*0.5 0+(dv-1)*0.05 0+(dv-1)*0.05])
                hold on
            end

        end
        
        ylabel("Zscored Frequency",'FontSize',15); 
        xlabel({"Speed (cm/s)"},'FontSize',20); 
        title(strcat(samplename,", All cell"),'FontSize',20);
        xticklabels(["0-5" "5-10" "10-15" "15-20"])
        ax = gca;
        ax.XAxis.FontSize = 12
        ax.YAxis.FontSize = 12
        exportgraphics(gcf,strcat(CDir,"\ZScored_all_trial\ZScored_DVspeed_",samplename,".jpg"))
clf
    end
end

save("ZSpeedHz.mat","gp2_ZSpeedHz","gp2_SpeedHz",'-append')
 
