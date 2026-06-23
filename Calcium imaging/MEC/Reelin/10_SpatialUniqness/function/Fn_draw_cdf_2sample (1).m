function Fn_draw_cdf_2sample(DATA1,DATA2,XLABEL,SAMPLELABELS,TITLE)


% DATA1=wt
% DATA2=eb


h1=cdfplot(DATA1);
hold on
h2=cdfplot(DATA2);
h1.LineWidth=3;h2.LineWidth=3;
h1.Color="black"
xlabel(XLABEL,'FontSize',15);
ylabel("Proportion",'FontSize',15);
title(TITLE);
yticks([0 0.5 1]);
box off

hold on
[f1,x1] = ecdf(DATA1);
[f2,x2] = ecdf(DATA2);
Idx = knnsearch(x2,x1);
ksd=abs(f2(Idx)-f1);
ksd_id=find(ksd==max(ksd));
plot([x1(ksd_id) x2(Idx(ksd_id))],[f1(ksd_id) f2(Idx(ksd_id))],'k','LineWidth',2);

grid off
legend(SAMPLELABELS(1),SAMPLELABELS(2),'Location','southeast','FontSize',15);
% title({TITLE;strcat("p = ",num2str(p), ", KS distance = ", num2str(ks2stat) )})


[h,p,ks2stat] = kstest2(DATA1, DATA2);
if p>0.0001
    title({TITLE;strcat("KS test, p = ", num2str(round(p,4)), ...
        ", KS distance = ", num2str(round(ks2stat,4)))})
else
    title({TITLE;strcat("KS test, p<0.0001, KS distance = ", num2str(round(ks2stat,4)) ) })  
end
    
end