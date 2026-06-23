function Fn_draw_histogram_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE,BIN)


% DATA1=wt
% DATA2=eb


h1 = histogram(DATA1,"BinWidth",BIN,"Normalization","probability");
hold on
h2 = histogram(DATA2,"BinWidth",BIN,"Normalization","probability");
h1.LineWidth=2;
h2.LineWidth=2;
h1.FaceColor="black";
h1.EdgeColor="black";
xlabel(XLABEL,'FontSize',15);
ylabel("Probability",'FontSize',15);
title(TITLE);
% yticks([0 0.5 1]);


plotLabels = string(SAMPLELABELS);
legend(plotLabels(1),plotLabels(2),'Location','northwest','FontSize',15);
[~,p,~] = kstest2(DATA1, DATA2);
if p>0.0001
    title({TITLE;strcat("KS test, p = ", num2str(round(p,4)) )})
else
    title({TITLE;strcat("KS test, p<0.0001" )})  
end


    
end
