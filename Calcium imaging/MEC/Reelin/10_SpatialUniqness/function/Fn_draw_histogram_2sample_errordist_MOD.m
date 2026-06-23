function Fn_draw_histogram_2sample_errordist_MOD(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE)


% DATA1=wt
% DATA2=eb


histogram(DATA1,"BinWidth",5,"Normalization","probability");
hold on
histogram(DATA2,"BinWidth",5,"Normalization","probability");
h1.LineWidth=2;
h2.LineWidth=2;
h1.Color="black";
xlabel(XLABEL,'FontSize',15);
ylabel("Probability",'FontSize',15);
title(TITLE);
% yticks([0 0.5 1]);
ylim([0 0.1]);

legend(SAMPLELABELS(1),SAMPLELABELS(2),'Location','northeast','FontSize',15);
[h,p,ks2stat] = kstest2(DATA1, DATA2);
if p>0.0001
    title({TITLE;strcat("KS test, p = ", num2str(round(p,4)) )});
else
    title({TITLE;strcat("KS test, p<0.0001" )})  ;
end


    
end