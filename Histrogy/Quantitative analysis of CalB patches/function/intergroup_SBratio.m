function I_SBratio = intergroup_SBratio(dataS, dataBG, Snum, figxlim, figylim, filenameS)
%% figure handling number
fig = gcf;
i_fig = fig.Number;
if i_fig > 1; i_fig = i_fig+1; end

 meanback = cell(Snum,1);

%% background division
  
  I_SBratio = cell(Snum,1)

for j = 1:3:Snum

    for ii=1:1:3
        i = ii+j-1;
     
        meanback{i} = movmean(dataBG{i}(:,2),100)
        I_SBratio{i} = dataS{i}(:,2) ./  meanback{i} 


        subplot(3,1,ii)
        plot(dataS{i}(:,1), I_SBratio{i});
        hold on;
        xlim(figxlim);
        ylim(figylim); % Change y lim
        title(strcat("signal ratio (intensity / mean intensity), ", erase(filenameS{i}, ".csv")));

    end
    xlabel('distance [um]');
    ylabel('ratio');
    subplot(3,1,1)
    print(strcat("results/SBratio ", erase(filenameS{j}, ".csv"), ".jpg"), '-djpeg', '-r300')
    close all
end
disp('div fin')

end
