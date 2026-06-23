function theta = calculate_angle(A, B)

    delta_y = A(:,2)*-1 - B(:,2)*-1;
    delta_x = A(:,1) - B(:,1);
    

    theta_rad = atan2(delta_y, delta_x);
    

    theta = rad2deg(theta_rad);
end
