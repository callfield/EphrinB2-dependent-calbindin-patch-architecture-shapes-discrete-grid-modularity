
function [AngleArray, newV] = extractAngle(DLC_data, RecStart, newV, x1, y1, x2, y2, AngleArray)

    A = DLC_data(:, [x1, y1]);
    B = DLC_data(:, [x2, y2]);
    theta = calculate_angle(A, B);
    load("ST_dF_grid_aut_data.mat", "sTrk");

    for k = 1:size(RecStart, 1)
        tmp=theta(RecStart(k,1):RecStart(k,2));
        % Remove frames with invalid tracking coordinates.
        susp=find(sTrk{newV}(:,2)==0&sTrk{newV}(:,3)==0&sTrk{newV}(:,4)==0);%
        tmp(susp,:)=[];
        AngleArray = [AngleArray; tmp];
        newV=newV+1;
    end



    
end
