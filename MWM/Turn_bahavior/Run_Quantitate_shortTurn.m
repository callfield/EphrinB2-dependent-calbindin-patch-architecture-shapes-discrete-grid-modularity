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
% Group1=[1:2:20]; % animal ID for group1
% Group2=[2:2:20]; % animal ID for group1

%%




Angle=cell(numAnimal);
AngVelocity=cell(numAnimal);
for i=1:numAnimal

    x=movmean(Track{i}(:,3),5);
    y=movmean(Track{i}(:,4),5);
    refpos=[0,0,0];
    pos=[x y];
    pos(:,3)=0;
    pos=pos.';
    Angle{i}=zeros(length(x),1);
    for ff=2:length(x)
        [tmp_rng,tmp_ang] = rangeangle(pos(:,ff)-pos(:,ff-1));
        Angle{i}(ff)=tmp_ang(1);

    end
    tmp=diff(Angle{i});
    tmp(find(tmp>90))=tmp(find(tmp>90))-360;
    tmp(find(tmp<-90))=tmp(find(tmp<-90))+360;
    AngVelocity{i}=tmp;


end


ST=4;
EN=25*59; % 25 Hz



% find crossing position cos=0


for i=1:numAnimal
    
    xy1=[];
    xy1(1,:)=ST:EN;
    xy1(2,:)=cos(deg2rad(Angle{i}(ST:EN)));

    XY2=[];
    XY2(1,:)=ST:EN;
    XY2(2,:)=cos(deg2rad(Angle{i}(ST)));
    if cos(deg2rad(Angle{i}(ST)))<-0.95 || cos(deg2rad(Angle{i}(ST)))>0.95
        XY2(2,:)=cos(deg2rad(Angle{i}(ST)))*0.95;
    end
    [Cross_frame{i},y,uc_code]=intersection_of_two_curves(xy1(1,:),xy1(2,:),XY2(1,:),XY2(2,:));
end


%% Define turn based on crossing of trace



WIN=5; % moving window

Cross_turn_frame=cell(numAnimal,1);
for i=1:numAnimal


    x=Track{i}(:,3);
    y=Track{i}(:,4);

    ST=1;
    END=ST+WIN;
    ST2=END+1;
    END2=ST2+WIN;

    TURN=1;
    while END < (length(x)-30)
        
        xy=Track{i}(ST:END,3:4).';
        XY=Track{i}(ST2:END2,3:4).';
       [crossX,crossY,uc_code]=intersection_of_two_curves(xy(1,:),xy(2,:),XY(1,:),XY(2,:));


       if isempty(crossX)==1
           END=END+WIN;
       elseif isempty(crossX)==0 & length(crossX)==1
          
           CrossPoint=[crossX crossY];
           cross_ST=min(knnsearch(xy.',CrossPoint,'K',1));
           cross_ED=min(knnsearch(XY.',CrossPoint,'K',1));

           Cross_turn_frame{i}(TURN,1)=cross_ST+ST-1;
           Cross_turn_frame{i}(TURN,2)=cross_ED+ST2-1;

           ST=cross_ED+ST2;
           END=ST+WIN;

           TURN=TURN+1;
       elseif isempty(crossX)==0 & length(crossX)>1
           i
           crossX
           ST
           END=END+WIN
       end
      
        ST2=END+1;
        END2=ST2+WIN;

    end


end


% Distance during turn
for i=1:numAnimal
   for t=1:size(Cross_turn_frame{i},1)
        TURN_st=Cross_turn_frame{i}(t,1);
        TURN_ed=Cross_turn_frame{i}(t,2);

        x=Track{i}(TURN_st:TURN_ed,3);
        y=Track{i}(TURN_st:TURN_ed,4);


%  distance from center of turn to center of arena 
        tX=median(x,"all","omitnan");
        tY=median(y,"all","omitnan");
        MedDistance_Center=sqrt(tX^2+tY^2);

% average distance during turn on from plat form
        AveDistance_target=mean(Track{i}(TURN_st:TURN_ed,24),"all","omitnan");
        MinDistance_target=min(Track{i}(TURN_st:TURN_ed,24),[],"omitnan");

% Swiming distance on each turn
        x=Track{i}(TURN_st:TURN_ed,3);
        y=Track{i}(TURN_st:TURN_ed,4);
        dX=diff(x);
        dY=diff(y);
        SwimmDistance=sum(sqrt(dX.^2+dY.^2));

        Cross_turn_frame{i}(t,3)=MedDistance_Center;
        Cross_turn_frame{i}(t,4)=AveDistance_target;
        Cross_turn_frame{i}(t,5)=MinDistance_target;
        Cross_turn_frame{i}(t,6)=SwimmDistance;

   end
end




%%

OmitLongTurn=Cross_turn_frame;
for i = 1:numAnimal
% omit longer than 3/4 of circumference and short turn
        OmitLongTurn{i}(find(OmitLongTurn{i}(:,6)>150*pi*3/4),:)=[];
        OmitLongTurn{i}(find(OmitLongTurn{i}(:,6)<5),:)=[];



end





%% visualise

COLOR=[0.9 0.9 0 ; 0.9 0 0 ; 0 0.9 0 ; 0 0 0.9];

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
    h = plot(xunit, yunit,"LineWidth",2,"Color",[.6 .6 .6]);


if size(OmitLongTurn{i},1)<4
   turnNum=size(OmitLongTurn{i},1);
else
   turnNum=4;
end


    x=Track{i}(:,3);
    y=Track{i}(:,4);
hold on
plot(x(1 :OmitLongTurn{i}(turnNum,2) ),y(1 :OmitLongTurn{i}(turnNum,2)),'k--',"LineWidth",1.5)
for t=1:turnNum
plot(x(OmitLongTurn{i}(t,1) :OmitLongTurn{i}(t,2) ),y(OmitLongTurn{i}(t,1) :OmitLongTurn{i}(t,2)),"LineWidth",3,"Color",COLOR(t,:))
end
    title(strcat("Animal",num2str(AnimalID(i))))
    xlim([-r r]);
    ylim([-r r]);
    pbaspect([1 1 1]);
    exportgraphics(gcf,strcat("OmitLongTurn crossTurn4 Probe Track animal",num2str(AnimalID(i)) ,".jpg"))
end
clf


%%

% Distance on each turn
MinDistance_ToPlatform=nan(22,4);
for a = 1:22

        for t=1:size(OmitLongTurn{a}(:,4),1)
            if t<5
             MinDistance_ToPlatform(a,t)=OmitLongTurn{a}(t,5);
            end
        end   
    
end

writematrix(MinDistance_ToPlatform(Group1,:),"MinDistance_ToPlatform.xlsx","Sheet","Group1")
writematrix(MinDistance_ToPlatform(Group2,:),"MinDistance_ToPlatform.xlsx","Sheet","Group2")




%%

% In zone ratio
InZoneRatio=cell(4,1);
for i = 1:numAnimal
    if size(OmitLongTurn{i},1)<3
       turnNum=size(OmitLongTurn{i},1);
    else
       turnNum=3;
    end


    for t=1:turnNum
        TURN_st=OmitLongTurn{i}(t,1);
        TURN_ed=OmitLongTurn{i}(t,2);
        TAR=Track{i}(TURN_st:TURN_ed,19);% In zone target quad
        Q1=Track{i}(TURN_st:TURN_ed,21);% In zone 1st quad
        Q2=Track{i}(TURN_st:TURN_ed,22);% In zone 2nd quad
        Q3=Track{i}(TURN_st:TURN_ed,23);% In zone 3rd quad

        InZoneRatio{1}(i,t)=length(find(TAR==1))/length(TAR);
        InZoneRatio{2}(i,t)=length(find(Q1==1))/length(Q1);
        InZoneRatio{3}(i,t)=length(find(Q2==1))/length(Q2);
        InZoneRatio{4}(i,t)=length(find(Q3==1))/length(Q3);

    end
end


writematrix(InZoneRatio{1}(Group1,:),"Ratio_Target.xlsx","Sheet","Group1")
writematrix(InZoneRatio{1}(Group2,:),"Ratio_Target.xlsx","Sheet","Group2")



%Cumulative In zone ratio
CumSum_InZoneRatio=cell(4,1);
for i = 1:numAnimal
    if size(OmitLongTurn{i},1)<3
       turnNum=size(OmitLongTurn{i},1);
    else
       turnNum=3;
    end

TAR=[];
Q1=[];
Q2=[];
Q3=[];

    for t=1:turnNum
        TURN_st=OmitLongTurn{i}(t,1);
        TURN_ed=OmitLongTurn{i}(t,2);

        tmp=Track{i}(TURN_st:TURN_ed,19);% In zone target quad
        TAR=[TAR;tmp];

        tmp=Track{i}(TURN_st:TURN_ed,21);% In zone 1st quad
        Q1=[Q1;tmp];

        tmp=Track{i}(TURN_st:TURN_ed,22);% In zone 2nd quad
        Q2=[Q2;tmp];

        tmp=Track{i}(TURN_st:TURN_ed,23);% In zone 3rd quad
        Q3=[Q3;tmp];

       CumSum_InZoneRatio{1}(i,t)=length(find(TAR==1))/length(TAR);
       CumSum_InZoneRatio{2}(i,t)=length(find(Q1==1))/length(Q1);
       CumSum_InZoneRatio{3}(i,t)=length(find(Q2==1))/length(Q2);
       CumSum_InZoneRatio{4}(i,t)=length(find(Q3==1))/length(Q3);
     end


    
end

writematrix(CumSum_InZoneRatio{1}(Group1,:),"CumsSumRatio_Target.xlsx","Sheet","Group1")
writematrix(CumSum_InZoneRatio{1}(Group2,:),"CumsSumRatio_Target.xlsx","Sheet","Group2")








