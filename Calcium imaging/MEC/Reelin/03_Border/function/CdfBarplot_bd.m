function CdfBarplot_bd(OutDir,data1,data2,name1, name2, xlbl,ylbl,ttl,filename)
clf
h1=cdfplot(data1)
hold on
h2=cdfplot(data2)
set(gca, 'FontSize',10)
h1.LineWidth=2;h2.LineWidth=2
h1.Color="black"
xlabel(xlbl,'FontSize',18)
ylabel(ylbl,'FontSize',20);
yticks([0 0.5 1]);
box off

hold on
[f1,x1] = ecdf(data1);
[f2,x2] = ecdf(data2);
Idx = knnsearch(x2,x1);
ksd=abs(f2(Idx)-f1);
ksd_id=find(ksd==max(ksd));
plot([x1(ksd_id) x2(Idx(ksd_id))],[f1(ksd_id) f2(Idx(ksd_id))],':k','LineWidth',2)

grid off
legend(name1,name2,'Location','southeast','FontSize',15)
[h,p,ks2stat] = kstest2(data1, data2)
title({ttl;strcat("p = ",num2str(p),...
    ", KS distance = ", num2str(round(ks2stat,2) ),...
    ", ",name1,"=", num2str(length(data1)),", ",name2,"=",num2str(length(data2)))})
exportgraphics(gcf,strcat(OutDir,"CDF ", filename),'Resolution',300)



% bar plot
clf
data1(isnan(data1))=[];
data2(isnan(data2))=[];
er1=std(data1)/length(data1)^0.5;
er2=std(data2)/length(data2)^0.5;

b=bar([mean(data1) mean(data2)],'FaceColor','flat','LineWidth',1)
b.CData(1,:) = [1 1 1];
b.CData(2,:) = [1 0.7 0.7];
hold on
er = errorbar(1:2,[mean(data1,'omitnan') mean(data2,'omitnan')],[er1 er2],'LineWidth',2,'CapSize',20);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
% ylim([0 max(er.YData)*1.3])
xticks([1 2]);
xticklabels([name1 name2])
%  set(gca, 'FontSize',10)
ax = gca;
ax.XAxis.FontSize = 16
[h,p]=ttest2(data1, data2);
 xlim([0.25 2.75]);
    ylabel(xlbl,'FontSize',15);
    title({strcat(ttl); ...
        strcat("ttest p= ", num2str(p))})
   exportgraphics(gcf,strcat(OutDir,"Bar ", filename),'Resolution',300)
close

end
