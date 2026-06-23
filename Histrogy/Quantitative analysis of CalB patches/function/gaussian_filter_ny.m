function I_gaussianfilt = gaussian_filter_ny(I_SBratio, winsize, Snum, figxlim, figylim, dataS, filenameS)
fig = gcf;
i_fig = fig.Number;
if i_fig > 1; i_fig = i_fig+1; end

figure(i_fig);

%% gaussian filter

% filter coefficients
b = gausswin(winsize);
a = sum(b);

I_gaussianfilt = cell(Snum,1);


for j = 1:3:Snum

    for ii=1:1:3
        i = ii+j-1;
        subplot(3,1,ii)
        x = dataS{i}(:,1);
        I_gaussianfilt{i} = filter(b,a, I_SBratio{i});

        plot(x, I_gaussianfilt{i}); hold on;

        xlim(figxlim);
        ylim(figylim);
        title(strcat("signal ratio (intensity / mean intensity) after gaussian filter, ", erase(filenameS{i}, ".csv")));
                                 % Correct Amplitudes Of Peaks
    end
    xlabel('distance [um]');
    ylabel('normalized intensity');
    subplot(3,1,1)
    print(strcat("results/Gfiltered_SBratio ", erase(filenameS{j}, ".csv"), ".jpg"), '-djpeg', '-r300')
    close all
end



end