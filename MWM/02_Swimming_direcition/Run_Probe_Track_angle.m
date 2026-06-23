clear all


%%% Required input
% Track{i} :  track of animal i (processed from ExoVision output csv)
% 1: Trial time [s], 2: Recording time [s], 3: X center [cm], 4: Y center [cm], 
% 5: X nose [cm], 6: Y nose [cm], 7: X tail [cm], 8: Y tail [cm], 
% 9: Area [cm²], 10: Area change [cm²], 11: Elongation, 12: Direction [deg], 
% 13: Latency to platform (Area / Center-point), 14: Latency to platform (Area / Nose-point), 15: Latency to platform (Area / Tail-base), 
% 16: Velocity (Center-point) [cm/s], 17: Velocity (Nose-point) [cm/s], 18: Velocity (Tail-base) [cm/s], 
% 19: In zone (Target Quad / Center-point), 20: In zone (Target Platform / Center-point), 21: In zone (Quad 1 / Center-point), 22: In zone (Quad 2 / Center-point), 23: In zone (Quad 3 / Center-point), 
% 24: Distance to point [cm], 25: Heading to point (Target Platform [Center] / Nose-point) [deg], 26: Heading to point (Target Platform [Center] / Tail-base) [deg], 27: Turn angle [deg], 
% 28: Angular velocity [deg/s]

% numAnimal=length(1:20)


%%

p=cell(numAnimal,1);
for i=1:numAnimal
    [time col]=size(Track{i});
    clf
    % box off
    axis off
    
    hold on
 r=75;
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) ;
    yunit = r * sin(th) ;
    h = plot(xunit, yunit,"LineWidth",5,"Color",[.6 .6 .6]);


   st=1;
   ed=round(time*5/60);
    x=Track{i}(st:ed,3);
    y=Track{i}(st:ed,4);

p{i} = polyfit(x,y,1);
y1 = polyval(p{i},x);
plot(x,y1,"LineWidth",2,"Color","red")

    plot(x,y,"LineWidth",3,"Color","black")
    
    title(strcat("Animal",num2str(AnimalID(i))))
    xlim([-r r]);
    ylim([-r r]);
    pbaspect([1 1 1]);


end
clf

Angle=nan(1,1);
for i=1:numAnimal
    Angle(i,1)=180*atan2(1,p{i}(1))/pi
end
writematrix(Angle,'Angel_first5.xlsx')
