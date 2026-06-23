%% http://nmarkou.blogspot.com/2012/03/cubic-spline-interpolation.html

function [ yy,xx ] = cubicSpline2d(N, points, gradients )
% CUBICSPLINE - returns N interpolation points using cubic spline
%   interpolation method, this method cann be used for space curves
%
% Input
% N         : number of interpolation points 
% points    : given points [x y]
% gradients : gradient at first and last point 
%
% Output    
% xx        : uniform spaced function interpolation at N points
% yy        : uniform spaced N points
%
%-------------------------------------------------------------------------
    %% Validate input arguments
    if isempty(N) || N < 1
        error('N must be >= 1');
    end
    
    if isempty(points) || size(points,1) < 1
        error('point must have >= 1 rows');
    end
    
    if isempty(points) || size(points,2) ~= 2
        error('point must have 2 collumns');
    end
    
    if isempty(gradients) || numel(gradients) ~= 4
        error('gradients must have 4 elements');
    end
    %%
    % get number of points
    [rows ~] = size(points);
    
    % get total length of points
    u = [1 : rows]';
    x = [points(:,1)];
    y = [points(:,2)];
    
    [xx,~] = cubicSpline(N, [u x],gradients(:,1));
    [yy,~] = cubicSpline(N, [u y],gradients(:,2));
end