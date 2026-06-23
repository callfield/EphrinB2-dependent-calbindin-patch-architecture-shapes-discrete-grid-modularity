function Fn_draw_bar_2sample(DATA1,DATA2,XLABELS,YLABEL,TITLE)


% DATA1=wt
% DATA2=eb

DATA1(isnan(DATA1))=[];
DATA2(isnan(DATA2))=[];
er1=std(DATA1)/length(DATA1)^0.5;
er2=std(DATA2)/length(DATA2)^0.5;

b=bar([mean(DATA1) mean(DATA2)],'FaceColor','flat','LineWidth',1)
b.CData(1,:) = [1 1 1];
b.CData(2,:) = [1 0.7 0.7]
hold on
% scatter([repmat(1,length(data1),1);repmat(2,length(data2),1)],[data1;data2],'k',"filled")
er = errorbar(1:2,[mean(DATA1,'omitnan') mean(DATA2,'omitnan')], ...
    [er1 er2],'LineWidth',2,'CapSize',20);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
% ylim([1 1.8])
xticks([1 2]);
xticklabels(XLABELS)

ax = gca;
ax.XAxis.FontSize = 16
[h,p]=ttest2(DATA1, DATA2);
 xlim([0.25 2.75]);%ylim([50 300])
    ylabel(YLABEL,"FontSize",17);
    title(TITLE,"FontSize",15)
text(1,1.7,strcat("ttest p= ", num2str(round(p,4))),"FontSize",12)
    
end