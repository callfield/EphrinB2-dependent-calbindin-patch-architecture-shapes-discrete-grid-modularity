function [ yy, xx] = cubicSpline(N,points,gradients)
% CUBICSPLINE - returns N interpolation points using cubic spline
%   interpolation method, this method cannot be used for space curves
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
    
    if isempty(gradients) || numel(gradients) ~= 2
        error('gradients must have 2 elements');
    end
    %% coefficient calculation part
    
    % get number of points
    [rows ~] = size(points);
    
    % compute inverse matrix to be used
    matrix = inv([1 0 0 0 ; 0 1 0 0 ; 1 1 1 1; 0 1 2 3]);
    
    % initialize coefficients structure
    coefficients = zeros(rows-1,4);
    
    % given n points we must calculate n-1 polynomials
    for i = 2 : rows
        pEnd = [];
        pStart = [];
        
        % calculate gradient using finite central differences
        if (i-1) == 1
            pStart = gradients(1);
        else
            pStart = (points(i,2) - points(i-2,2))/2;
        end
        
        if i == rows
            pEnd = gradients(2);
        else
            pEnd = (points(i+1,2) - points(i-1,2))/2;
        end
        
        % create vector [Pi P'i Pi+1 P'i+1]'
        vector = [points(i-1,2);pStart;points(i,2);pEnd];
        
        % calculate polynomial coefficients
        coefficients(i-1,:) = (matrix * vector)';
    end
    
    %% interpolation part
    
    % get max X and min X and interval
    minX = points(1,1);
    maxX = points(end,1);
    intervalX = (maxX - minX) / (N - 1);
    xx = minX : intervalX : maxX;
    
    % interpolate at given locations
    yy = zeros(1,N);
    splineIndex = 1;
    for i = 2 : N-1
        x = xx(i);
        % find the index of the used spline
        for j = splineIndex : rows
            if x >= points(j,1) && x < points(j+1,1)
                splineIndex = j;
                break;
            end
        end
        splineCoeffs = coefficients(splineIndex,:);
        % compute m 
        m = (xx(i) - points(splineIndex,1))/...
            (points(splineIndex+1,1) - points(splineIndex,1));
        % compute value with given spline and m
        yy(i) = splineCoeffs(1) + splineCoeffs(2) * m + ...
            splineCoeffs(3) * m^2 + splineCoeffs(4) * m^3;
    end
    yy(1) = points(1,2);
    yy(end) = points(end,2);
end
