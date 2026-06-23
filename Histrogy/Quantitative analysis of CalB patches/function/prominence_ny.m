function [pks,locs,w,p_from_varey, wxPk] =  prominence_ny(I_gaussianfilt, Snum, figxlim, figylim, dataS, cutlength, minpeakprominence,filenameS)
% pks: peak intensity
% locs: peak location in x axis
% w: width of half maximum
% p_from_varey: peak intensity from base
% wxPk: location of half-width maximum


fig = gcf;
i_fig = fig.Number;
if i_fig > 1; i_fig = i_fig+1; end

figure(i_fig);


%% prominence


pks = cell(Snum, 1); locs = cell(Snum, 1); w = cell(Snum, 1); p_from_varey = cell(Snum, 1);
for j = 1:3:Snum

    for ii=1:1:3
        i = ii+j-1;
        subplot(3,1,ii)
        x = dataS{i}(cutlength:end,1);
        y = I_gaussianfilt{i}(cutlength:end,1);

        plot(x, y); hold on;
        [pks{i},locs{i},w{i},p_from_varey{i},wxPk{i}] =  findpeaks_ho(y,x,'Annotate','extents', 'MinPeakProminence', minpeakprominence);% save variables for prominence analysis

        findpeaks(y,x,'Annotate','extents', 'MinPeakProminence',minpeakprominence); % plot prominence analysis
        text(locs{i}-100, pks{i}+1,num2str(p_from_varey{i}));
        text(locs{i}-100, pks{i}-1,num2str(w{i}));

        ylim(figylim);
       % ylim([0 5]);

        legend('hide');
        title(strcat("peak analysis, ", erase(filenameS{i}, ".csv")));

    end
    xlabel('distance [um]');
    ylabel('intensity ratio');
    subplot(3,1,1)
    print(strcat("results/peak ", erase(filenameS{j}, ".csv"), ".jpg"), '-djpeg', '-r300')
    close all
end



% end