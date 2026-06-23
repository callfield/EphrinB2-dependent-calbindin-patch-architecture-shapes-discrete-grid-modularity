   close all; clear all

addpath('function');
addpath(pwd);

load('../08_GridModule/GridMod.mat'); 

load('..\Data.mat'); 

load('ZSpeedHz_each1.mat');

load('..\Data.mat');

CDir=pwd;


%%
gp1_ZSpeedHz_fit_mean=cell(2,1);
ii=1;
for s=1:5
    for t=1:3
        if isnan(gp1_GSratio(s,t)) ==0


            gcell= find(GROUP1{s,t}(:,4)~=0);
             tmp2=cell(2,1);
            for m=1:2
                tmp_ModGcell=gcell(gp1_Mod{s,t}{m});
    
                tmp=gp1_ZSpeedHz{s,t}(1:20,tmp_ModGcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean
                hold on

                for l=1:size(tmp,2)
                    x=1:size(tmp,1);
                    y=tmp(1:size(tmp,1),l);
                    [p,S] = polyfit(x,y,1); 
                    [y_fit,delta] = polyval(p,x,S);
                    if m==1
                        plot(x,y_fit,'k-')  
                    else
                        plot(x,y_fit,'r-')
                    end
                    tmp2{m}(l)=p(1);
                end
              
                 gp1_ZSpeedHz_fit_mean{m}(ii)=nanmean(tmp2{m},2);
            
            end
            ii=ii+1;
        end
    end
end

gp1_ZSpeedHz_fit_mean{1}
gp1_ZSpeedHz_fit_mean{2}

gp2_ZSpeedHz_fit_mean=cell(2,1);
ii=1;
for s=1:7
    for t=1:3
        if isnan(gp2_GSratio(s,t)) ==0

            gcell= find(GROUP2{s,t}(:,4)~=0);
             tmp2=cell(2,1);
            for m=1:2
                tmp_ModGcell=gcell(gp2_Mod{s,t}{m});
    
                tmp=gp2_ZSpeedHz{s,t}(1:20,tmp_ModGcell);
                tmp(tmp==0)=NaN;% convert zero to nan before mean
                hold on
                for l=1:size(tmp,2)
                    x=1:size(tmp,1);
                    y=tmp(1:size(tmp,1),l);
                    [p,S] = polyfit(x,y,1); 
                    [y_fit,delta] = polyval(p,x,S);
                    if m==1
                        plot(x,y_fit,'k-')  
                    else
                        plot(x,y_fit,'r-')
                    end
                    
                    tmp2{m}(l)=p(1);
                end
              
                 gp2_ZSpeedHz_fit_mean{m}(ii)=nanmean(tmp2{m},2);
            
            end
            ii=ii+1;
        end
    end
end

Dif_gp1_ZSpeedHz_fit_mean=gp1_ZSpeedHz_fit_mean{1}-gp1_ZSpeedHz_fit_mean{2};
Dif_gp2_ZSpeedHz_fit_mean=gp2_ZSpeedHz_fit_mean{1}-gp2_ZSpeedHz_fit_mean{2};

writecell(gp1_ZSpeedHz_fit_mean, 'SpeedLine_fit.xlsx','Sheet', 'gp1_mod',"AutoFitWidth",false)
writecell(gp2_ZSpeedHz_fit_mean, 'SpeedLine_fit.xlsx','Sheet', 'gp2_mod',"AutoFitWidth",false)


writematrix(Dif_gp1_ZSpeedHz_fit_mean, 'SpeedLine_fit.xlsx','Sheet', 'Dif_gp1_mod',"AutoFitWidth",false)
writematrix(Dif_gp2_ZSpeedHz_fit_mean, 'SpeedLine_fit.xlsx','Sheet', 'Dif_gp2_mod',"AutoFitWidth",false)





