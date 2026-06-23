
% intersection_of_two_curves.m
%
% User-defined function that returns the intersection coordinates of two
% planar curves represented as polylines.
%
% Verified with MATLAB Home R2019a.
%
% Each curve may be a collection of multiple sub-curves joined into a single
% vector. In that case, insert NaN as a separator point between sub-curves.
% The function reports only intersections between curve 1 and curve 2.
% Intersections between sub-curves within the same curve are ignored.
%
% [x,y,uc_code]=intersection_of_two_curves(x1,y1,x2,y2)
%
% Inputs:
%   x1:      x-coordinates of the point sequence for curve 1 (row vector)
%   y1:      y-coordinates of the point sequence for curve 1 (row vector)
%            Usually, processing is faster if the curve with fewer vertices
%            is used as curve 1.
%   x2,y2:   x- and y-coordinates of curve 2, in the same format.
%
% Outputs:
%   x:       x-coordinates of intersections (row vector; one element per
%            intersection)
%   y:       y-coordinates of intersections
%   uc_code: Information for intersections whose detection was uncertain
%            (5 rows x number of uncertain cases).
%            Row 1: vertex index at the lower-numbered end of the uncertain
%                   segment in curve 1.
%            Row 2: vertex index at the lower-numbered end of the uncertain
%                   segment in curve 2.
%            Rows 3 and 4: which end of each segment caused the uncertainty:
%                   0: neither end,
%                   1: lower-numbered end,
%                   2: higher-numbered end,
%                   3: not an endpoint issue; the two segments overlap.
%            Row 5: intersection index if coordinates were output despite
%                   uncertainty; 0 if the intersection coordinates were not
%                   output.

function[x,y,uc_code]=intersection_of_two_curves(x1,y1,x2,y2)

% MATLAB did not appear to have a command for this exact task, so this
% function was written manually.
%
% Remaining issues: missed detections, duplicate detections, and ambiguous
% interpretations.
%
% If one curve passes exactly through a vertex of the other curve, an
% intersection may be missed or detected twice. If vertices from the two
% curves overlap, up to four duplicate detections may occur. If one segment
% from each curve overlaps perfectly, some users may prefer to interpret
% this as infinitely many intersections rather than no intersection.
%
% In the current implementation, ambiguous cases should be reviewed and
% corrected by the user based on uc_code.
%
% Remaining issue: processing time.
%
% Detection can be slow when a long segment in curve 1, especially one with
% a slope near 45 degrees, overlaps the search region of many short segments
% in curve 2. This is usually not a problem for curves with few segments, but
% larger datasets may need further optimization.

uc_code=[];   % Stores information for uncertain intersection detections.

if (size(x1,1)~=1 | size(y1,1)~=1 | size(y2,1)~=1 | size(y2,1)~=1)
  disp('Error: input data must be row vectors.')
  x=NaN;
  y=NaN;
  return
end

N1=length(x1);   % Length of the input curve vectors.
N2=length(x2);

if N1~=length(y1)||N2~=length(y2)
  disp('Error: x and y vectors must have the same length.')
  x=NaN;
  y=NaN;
  return
end

if ~isreal([x1 y1 x2 y2]);   % NaN is treated as real, so it is allowed.
  disp('Error: complex-valued data cannot be used.')
  x=NaN;
  y=NaN;
  return
end

ind2 = [1:N2];   % Assign indices to the vertices of curve 2 so original
                 % vertex numbers remain traceable after elements are removed.

% Relative tolerance reference for coordinate values. eps(mar) is used as
% the tolerance for judging vertex positions.
mar=max(abs([x1 y1 x2 y2]));
                 % Maximum absolute coordinate value among numeric elements
                 % of all input curve vertices. NaN elements are ignored.

x=[];            % Prepare output vectors for intersection coordinates.
y=[];

% Variables that record segment information when an intersection decision is
% uncertain ("UnCertain").
uc_n1=[];        % Lower-numbered vertex index of the uncertain segment in curve 1.
uc_n2=[];        % Lower-numbered vertex index of the uncertain segment in curve 2.
uc_x1=[];        % Which end of each segment caused the uncertainty.
uc_x2=[];        % 0: neither end, 1: lower-numbered end,
                 % 2: higher-numbered end, 3: overlapping segments.
uc_ot=[];        % Intersection index if an uncertain intersection was output;
                 % 0 if the intersection coordinates were not output.

for n1=[1:N1-1]; % Process each segment of curve 1 in order.

  % Extract both endpoints of the current segment in curve 1.
  % p and s stand for primary and secondary.
  % Note that the coordinate order of p and s is not guaranteed.
  x1p=x1(n1);    % x-coordinate of the lower-numbered endpoint.
  y1p=y1(n1);    % y-coordinate of the lower-numbered endpoint.
  x1s=x1(n1+1);  % x-coordinate of the higher-numbered endpoint.
  y1s=y1(n1+1);  % y-coordinate of the higher-numbered endpoint.
  if isnan(x1p+y1p+x1s+y1s)
    continue;    % If any coordinate is NaN, skip this segment. A segment
                 % containing NaN is not drawn and is treated as nonexistent.
  end

  x2t=x2;        % Copy curve 2 so it can be edited during processing.
  y2t=y2;
  ind2t=ind2;

  % If either the x or y coordinate of a curve-2 vertex is NaN, set both
  % coordinates of that vertex to NaN.
  idnan = ~( ~isnan(x2t) .* ~isnan(y2t) ); % 1 only for elements containing NaN.

  x2t(idnan==1) = NaN;
  y2t(idnan==1) = NaN;

  % To shorten processing time, keep only the curve-2 segments that need to
  % be checked.

  % Flags indicating whether each curve-2 polyline segment should be checked
  % for intersections. 1: required, 0: not required.
  need2 = ones(1,N2);        % Initially mark every segment as required.

  % Search region for the current segment of curve 1: an axis-aligned
  % rectangle containing the segment, with a small margin.
  st = max([y1p y1s]) + eps(mar);  % y-coordinate of the top edge.
  sb = min([y1p y1s]) - eps(mar);  % y-coordinate of the bottom edge.
  sl = min([x1p x1s]) - eps(mar);  % x-coordinate of the left edge.
  sr = max([x1p x1s]) + eps(mar);  % x-coordinate of the right edge.

  % Remove curve-2 segments that cannot intersect the search region. For
  % each of the four boundaries (top, bottom, left, right), check which side
  % of the boundary each curve-2 vertex lies on. Vertices on the opposite
  % side of the search region are tentatively marked as unnecessary.
  %
  % Even if a vertex is marked unnecessary, a neighboring required vertex may
  % form a segment crossing the search region. To avoid missing such
  % intersections, expand the required range by one vertex using expand_range.
  % Repeat this process for all four boundaries and finalize unnecessary
  % vertices.

  need2a = (y2t<st) & need2;  % Candidate removal: vertices above the region.
  need2a = expand_range(need2a);  % Rescue adjacent vertices (local function).
  need2b = (y2t>sb) & need2;  % Candidate removal: vertices below the region.
  need2b = expand_range(need2b);
  need2c = (x2t>sl) & need2;  % Candidate removal: vertices left of the region.
  need2c = expand_range(need2c);
  need2d = (x2t<sr) & need2;  % Candidate removal: vertices right of the region.
  need2d = expand_range(need2d);

  % Segments satisfying all four conditions are kept as candidates for
  % intersection search.
  need2 = need2a & need2b & need2c & need2d;

  % NaN separator vertices inside the search region may otherwise be marked
  % for deletion. If a NaN separator is removed, two intentionally separated
  % sub-curves may be interpreted as one connected curve, producing false
  % intersections on an artificial segment. To prevent this, keep all NaN
  % separator vertices.
  need2(isnan(x2t)==1) = 1;

  % Keeping every NaN separator can leave many unnecessary NaNs outside the
  % search region and increase processing time. Collapse consecutive NaN
  % vertices to a single NaN vertex.
  nan1=isnan(x2t);
  nan2=circshift(nan1,1);
  nan2(1)=0;
  need2((nan1 & nan2)==1) = 0;

  % Delete curve-2 vertices that are no longer needed.
  x2t(need2==0) = [];     % Remove x-coordinates for vertices with flag 0.
  y2t(need2==0) = [];     % Remove y-coordinates for vertices with flag 0.
  ind2t(need2==0) = [];   % Remove matching original vertex indices.
                          % This greatly reduces the number of segments to check.
  N2x=length(x2t);        % Number of curve-2 vertices after extraction.

  for n2=[1:N2x-1];       % Process only the extracted curve-2 segments.

    % Extract both endpoints of the current segment in curve 2.
    x2p=x2t(n2);
    y2p=y2t(n2);
    x2s=x2t(n2+1);
    y2s=y2t(n2+1);

    if isnan(x2p+y2p+x2s+y2s)
      continue;    % If any coordinate is NaN, skip this segment.
    end

    % Check the intersection of one segment from curve 1 and one segment
    % from curve 2.

    [xc,yc,pos_1,pos_2,uc_1,uc_2]= ...
                      line_cross(x1p,y1p,x1s,y1s,x2p,y2p,x2s,y2s,mar);

    % Note:
    % The naming of line identifiers differs between this function and the
    % local function. In this function, curve numbers are 1 and 2; in the
    % local function, the lines are a and b. Endpoint identifiers p and s in
    % this function correspond to endpoints 1 and 2 in the local function.

    cpout = 0;                % Flag indicating whether intersection coordinates
                              % have been output. Initially 0 (not output).

    if pos_1==0 && pos_2==0;  % Output coordinates only for true intersections.
                              % Intersections on extended lines are excluded.
      x=[x xc];
      y=[y yc];
      cpout = 1;              % Intersection coordinates have been output.
%     disp(['Detected intersections: ' num2str(length(x))]);     % Progress
%     disp(['  Segment position: ' num2str(n1) '/' num2str(N1)]); % display
%     lap=toc;                                                     % for the
%     disp(['  Elapsed time: ' num2str(lap) ' s']);               % command line.
    end

    % Record information when the intersection decision is uncertain.
    if ((uc_1~=0 | uc_2~=0) & ~isnan(uc_1));
      uc_n1=[uc_n1 n1];         % Curve-1 vertex index.
      uc_n2=[uc_n2 ind2t(n2)];  % Curve-2 vertex index.
        uc_x1=[uc_x1 uc_1];     % Which end of the curve-1 segment.
        uc_x2=[uc_x2 uc_2];     % Which end of the curve-2 segment
                                % (0: none, 1: start, 2: end, 3: overlap).
      if cpout==1;                 % If coordinates were output,
        uc_ot=[uc_ot length(x)];   % record the intersection element index.
      else
        uc_ot=[uc_ot 0];           % 0 if coordinates were not output.
      end
    end
  end;    % End scan over curve-2 segments.
end;    % End scan over curve-1 segments.

% Combine uncertainty information into one matrix for output.
uc_code=[uc_n1;uc_n2;uc_x1;uc_x2;uc_ot];

% =======================================================
% line_cross
%
% Local function that calculates the intersection of two lines.
%
% [x,y,pos_a,pos_b,uc_a,uc_b]= ...
%                   line_cross(xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2,xymax)
%
% Inputs:
%   xa1,ya1,xa2,ya2 : coordinates of both endpoints of segment a.
%   xb1,yb1,xb2,yb2 : coordinates of both endpoints of segment b.
%   xymax : maximum coordinate value, used as a reference for judging how
%           much segments a and b overlap.
%
% Outputs:
%   x,y : intersection coordinates.
%   pos_a,pos_b : where the intersection lies relative to segments a and b:
%                 0: on the segment,
%                 1: on the extension beyond endpoint 1,
%                 2: on the extension beyond endpoint 2.
%   uc_a,uc_b : whether the intersection decision is reliable for each
%               segment:
%                 0: reliable,
%                 1: uncertain near endpoint 1,
%                 2: uncertain near endpoint 2,
%                 3: segments a and b overlap.
%
% If the lines are parallel, or at least one segment has zero length (a
% point), x, y, pos_a, and pos_b are all NaN.

function [x,y,pos_a,pos_b,uc_a,uc_b] = ...
                     line_cross(xa1,ya1,xa2,ya2,xb1,yb1,xb2,yb2,xymax)

  % Let arbitrary points on segment a and segment b be:
  %     x = xa1+ka*(xa2-xa1)      x = xb1+kb*(xb2-xb1)
  %     y = ya1+ka*(ya2-ya1)      y = yb1+kb*(yb2-yb1)
  %     0 <= ka <= 1              0 <= kb <= 1
  %
  % Equating x and y on both lines gives:
  %     xa1 + ka*(xa2-xa1) = xb1 + kb*(xb2-xb1)
  %     ya1 + ka*(ya2-ya1) = yb1 + kb*(yb2-yb1)
  %     =>
  %     ka*(xa2-xa1) - kb*(xb2-xb1) = xb1-xa1
  %     ka*(ya2-ya1) - kb*(yb2-yb1) = yb1-ya1
  %     =>
  %     [xa2-xa1  xb1-xb2] [ka] = [xb1-xa1]
  %     [ya2-ya1  yb1-yb2] [kb]   [yb1-ya1]
  %
  % Compute ka and kb, then substitute them into the formulas above to get
  % the intersection coordinates. If ka or kb is outside the range above,
  % the true segments do not intersect; only their extended lines do.

  A=[xa2-xa1 xb1-xb2 ; ya2-ya1 yb1-yb2];
  B=[xb1-xa1;yb1-ya1];

  overlap = 0;   % Initial assumption: the two lines do not overlap.

  if rcond(A)<eps(1);  % The two segments are nearly parallel and unlikely
                       % to have an ordinary intersection.
                       % rcond(A) returns a small value when A is too
                       % ill-conditioned to invert accurately. It is a fairer
                       % criterion than det(A).

    % If the segments almost overlap, the decision is ambiguous and should
    % be flagged. First compute the distance between the two lines.
    aax=xa2-xa1;  % x component of segment-a vector.
    aay=ya2-ya1;  % y component of segment-a vector.

    % Let vector c run from endpoint 1 of segment a to endpoint 1 of segment b.
    ccx=xb1-xa1;  % x component of vector c.
    ccy=yb1-ya1;  % y component of vector c.

    % If a and b are parallel:
    %
    %                     abs(| a x c |)   | aax*ccy - aay*ccx |
    %     line spacing h = -------------- = ---------------------
    %                         | a |          sqrt(aax^2 + aay^2)
    %

    h = abs(aax*ccy - aay*ccx) / sqrt(aax^2 + aay^2 );   % Line spacing.

    x=NaN;               % Assume no intersection at first. The basic policy
    y=NaN;               % is to report no intersection for parallel lines,
    pos_a=NaN;           % regardless of the spacing, unless the overlap is
    pos_b=NaN;           % too close to ignore.
    if h<2*eps(xymax);   % If segments a and b almost overlap, forcing a
      overlap = 1;       % "no intersection" decision is questionable. Do
                         % not output coordinates, but flag the case for
                         % manual review.
    else           % If the lines are parallel and sufficiently separated,
      uc_a=NaN;    % no further checking is needed; return to the caller.
      uc_b=NaN;
      return;
    end
  end

  if overlap==1;   % Segments a and b almost overlap.
    uc_a=3;
    uc_b=3;
  else             % An intersection should exist; compute its coordinates.
    kk=inv(A)*B;   % Cases where an inverse cannot be formed were excluded,
    ka=kk(1);      % so this should not error here.
    kb=kk(2);
    x = xa1+ka*(xa2-xa1);   % Intersection coordinates.
    y = ya1+ka*(ya2-ya1);

    % If the intersection lies very close to a segment endpoint, mark the
    % decision as uncertain.
    % For segment a:
    uc_a=0;                 % Start with confidence in the decision.
    if abs(ka)<eps(1);      % Endpoint 1 is questionable.
      uc_a=1;               % Uncertain at endpoint 1.
    elseif abs(1-ka)<eps(1) % Endpoint 2 is questionable.
      uc_a=2;               % Uncertain at endpoint 2.
    end
    % For segment b: same treatment as for segment a.
    uc_b=0;
    if abs(kb)<eps(1)
      uc_b=1;
    elseif abs(1-kb)<eps(1)
      uc_b=2;
    end

    % Regardless of confidence, record where the crossing lies.
    if ka>=0 && ka<=1
      pos_a=0;      % Intersection lies on segment a.
    elseif ka<0
      pos_a=1;      % Intersection lies on the extension beyond endpoint 1.
    elseif ka>1
      pos_a=2;      % Intersection lies on the extension beyond endpoint 2.
    else
      pos_a=NaN;
    end
    if kb>=0 && kb<=1
      pos_b=0;      % Intersection lies on segment b.
    elseif kb<0
      pos_b=1;      % Intersection lies on the extension beyond endpoint 1.
    elseif kb>1
      pos_b=2;      % Intersection lies on the extension beyond endpoint 2.
    else
      pos_b=NaN;
    end
  end
end
% =======================================================

% =======================================================
% expand_range
%
% Local function that expands a processing target region by modifying a
% flag vector.
%
% [flg_new]=expand_range(flg1)
%
% Input:
%   flg1    : row vector of flag elements (0 or 1).
%
% Output:
%   flg_new : row vector containing the modified flags.
%             Modification: force any 0 element adjacent to a 1 element to 1.

function [flg_new]=expand_range(flg1)

  flg2 = circshift(flg1,-1);   % Shift the flag vector left.
  flg2(length(flg2)) = 0;      % Drop the left edge and insert 0 at the right.
  flg2 = flg1 | flg2;          % OR with the original vector, converting
                               % lower-numbered 0s adjacent to 1s into 1s.
  flg_new = circshift(flg2,1); % Then shift this vector right.
  flg_new(1) = 0;              % Drop the right edge and insert 0 at the left.
  flg_new = flg2 | flg_new;    % OR with the pre-shift vector, converting
                               % higher-numbered 0s adjacent to 1s into 1s.
end
% =======================================================

end        % End of intersection_of_two_curves.m

% ================================================================
% The following commented block is an evaluation script for this function.
% To run it, copy all lines below "% % xxtest_intersection_of_two_curves.m",
% remove the leading "% " from each line, save it as an appropriate .m file,
% and execute it. If needed, check MATLAB's character encoding with
% slCharacterEncoding().
% ================================================================

% % xxtest_intersection_of_two_curves.m
%
% % Test for the function intersection_of_two_curves.
%
% close all
% clear
%
% % =========================
% % Create model curve groups for evaluation (start).
% % Quadratic curves with various eccentricities.
%
% xd=-10;         % Directrix coordinate.
% xf=10;          % Focus coordinate.
% yy=[0:0.1:60];  % Equally spaced sample points on the y-axis.
% xyall=[];       % Stores the prototype of one curve group.
%
% % Scan over several eccentricities. Each curve becomes a sub-curve.
% for e=[0.2 0.38 0.52 0.62 0.71 0.78  0.88  1 ...
%                                  1.2 1.4 1.65 2 2.5 3.2 4.6 8 30]
%   if e==1;        % Parabola.
%     x = yy.^2/(2*(xf-xd)) + (xd+xf)/2;
%     xy = [x;yy];
%     xy2 = xy;
%     xy2(:,1) = [];
%     xy2 = fliplr(xy2);
%     xy2(2,:) = -xy2(2,:);
%     xy = [xy2 xy];
%   else
%     wx = e*(xf-xd)/(1-e^2);
%     wy = e*(xf-xd)/sqrt(abs(1-e^2));
%     bx = (xd*e^2-xf)/(1-e^2);
%     if e<1;       % Ellipse.
%       % Upper half of the curve.
%       x1 = -sqrt( wx^2*( 1 - (yy/wy).^2 )) - bx;  % Upper-left curve.
%       x2 =  sqrt( wx^2*( 1 - (yy/wy).^2 )) - bx;  % Upper-right curve.
%       xy = [ [x1 x2];[yy yy] ];          % Mixed left/right vertices.
%       nf = find( imag(xy(1,:)) ~= 0 );   % Points where the curve does not exist.
%       xy( :,nf ) = [];                   % Remove vertices with complex solutions.
%       xy = sortrows(xy');                % Sort mixed left/right vertices
%       xy = xy';                          % by ascending x-coordinate.
%       % Lower half of the curve.
%       xy2 = xy;                     % Copy the upper half.
%       xy2(:,1) = [];                % Remove the duplicated point,
%       xy2 = fliplr(xy2);            % sort by descending x-coordinate,
%       xy2(2,:) = -xy2(2,:);         % and invert the y-coordinate sign.
%       xy = [xy2 xy];                % Connect the upper and lower halves.
%     elseif e>1;   % Hyperbola.
%       % Left side.
%       x1 = -sqrt( wx^2*( 1 + (yy/wy).^2 )) - bx;  % Upper-left curve.
%       xy = [x1;yy];
%       xy2 = xy;                     % Copy this and create the lower half,
%       xy2(:,1) = [];                % as with the ellipse.
%       xy2 = fliplr(xy2);            % Connect the upper and lower curves.
%       xy2(2,:) = -xy2(2,:);         % Left side completed.
%       xy = [xy2 xy];
%       xym = xy;
%       % Right side.
%       x2 = sqrt( wx^2*( 1 + (yy/wy).^2 )) - bx;
%       xy = [x2;yy];
%       xy2 = xy;
%       xy2(:,1) = [];
%       xy2 = fliplr(xy2);
%       xy2(2,:)  =-xy2(2,:);         % Right side completed.
%       xy = [xy2 xy [NaN;NaN] xym];  % Left and right sides are separate
%                                     % curves, so connect them with NaN.
%     end
%   end
%   xyall = [xyall [NaN;NaN] xy];     % Connect sub-curves one after another
%                                     % through NaN vertices.
% end;      % Prototype for one curve group is complete.
%
% % Copy the prototype curve group and shift it to create curves 1 and 2.
% x1 = xyall(1,:) - 15;   % Curve 1.
% y1 = xyall(2,:) + 5;
% x2 = xyall(1,:) + 15;   % Curve 2.
% y2 = xyall(2,:) - 5;
%
% CC1 = [x1;y1];          % Matrix of x,y vertices for curve 1.
% CC2 = [x2;y2];          % Matrix of x,y vertices for curve 2.
% % Trim the outer boundary to remove unnecessary plotting range.
% CC1(:,(x1>80|x1<-80|y1>60|y1<-60)==1) = [];
% CC2(:,(x2>80|x2<-80|y2>60|y2<-60)==1) = [];
%
% % Open space for adding another curve.
% % Delete part of curve 1.
% CC1(:,(CC1(1,:)<-48 & CC1(2,:)>-30 & CC1(2,:)<25)==1) = NaN;
%
% % This creates many unnecessary NaN vertices, so collapse consecutive NaN
% % vertices to a single NaN vertex.
% nan1 = isnan(CC1(1,:));
% nan2 = circshift(nan1,1);
% nan2(1) = 0;
% CC1(:,(nan1 & nan2)==1) = [];
%
% % =============
% % Add another type of test curve: various double-line spacings and
% % parallelness values. On the plot, some differences are so small that the
% % lines appear to overlap.
%
% % Curve 1: twenty lines with slope 45 degrees.
% addo = [-75 -71;-28 -24];  % Line template.
% add1 = addo;               % Copy it and translate upward,
% Add1 = [];                 % placing 10 lines.
% for n=1:10
%   add1(2,:) = addo(2,:)+4*(n-1);
%   Add1 = [Add1 add1 [NaN;NaN]];
% end
% add1 = addo;               % Place another 10 lines to the right.
% add1(1,:) = addo(1,:) + 12;
% for n=1:10
%   add1(2,:) = addo(2,:) + 4*(n-1);
%   Add1 = [Add1 add1 [NaN;NaN]];
% end
% CC1 = [CC1 [NaN;NaN] Add1]; % Add these evaluation curves to curve group 1.
%
% % Curve 2: parallel lines with various spacings.
% add2 = addo;     % Draw each curve-2 line parallel to the 10 curve-1 lines
% Add2 = [];       % created on the left.
%                  % The lowest line overlaps curve 1 exactly.
%                  % The spacing increases slightly upward.
% for n=1:10
%   add2(2,:) = addo(2,:) + 4*(n-1) - eps(80)*(n-1);
%   Add2 = [Add2 add2 [NaN;NaN]];
% end
%
% % Curve 2: crossing lines with several tiny crossing angles.
% add2 = addo;     % Draw curve-2 lines through the centers of the 10
%                  % curve-1 lines created on the right, crossing at very
%                  % narrow angles.
%                  % The lowest line overlaps curve 1 exactly.
%                  % The crossing angle increases slightly upward.
% add2(1,:) = addo(1,:) + 12;
% for n=1:10
%   add2(2,1) = addo(2,1) + 4*(n-1) + eps(80)*(n-1)/45;
%   add2(2,2) = addo(2,2) + 4*(n-1) - eps(80)*(n-1)/45;
%   Add2 = [Add2 add2 [NaN;NaN]];
% end
% CC2 = [CC2 [NaN;NaN] Add2];  % Add these evaluation curves to curve group 2.
% % =============
%
% % Split curve coordinates into vectors for input to the test function.
% x1 = CC1(1,:);
% y1 = CC1(2,:);
% x2 = CC2(1,:);
% y2 = CC2(2,:);
%
% % Create model curve groups for evaluation (end).
% % =========================
%
% figure(1);          % Figure: many intersections between two curve groups.
% plot(x1,y1,'b',x2,y2,'r');                % Draw evaluation curve groups.
% axis equal
% hold on
% plot(x1,y1,'.b',x2,y2,'.r','MarkerSize',4);     % Add vertices to the plot.
%
% tic;                          % Start timer for processing time.
%
% % Call the user function being evaluated.
% [x,y,uc_code] = intersection_of_two_curves(x1,y1,x2,y2);
%
% delay=toc;                                      % Measure processing time.
%
% disp(' ')
% disp(' ')
% disp(' ')
% disp('[Uncertain information]')
% disp('  Segment number in curve 1')
% disp('  Segment number in curve 2')
% disp('  Problematic end of segment 1 (0: none, 1: start, 2: end')
% disp('  Problematic end of segment 2 (3: line overlap')
% disp('  Problematic intersection number; 0 if no intersection information')
% % Print the uncertainty code to the command line.
% uc_code
% disp(['Number of curve-1 vertices: ' num2str(length(x1))])
% disp(['Number of curve-2 vertices: ' num2str(length(x2))])
% disp(['Number of detected intersections: ' num2str(length(x))])
% disp(['Processing time: ' num2str(delay) ' s'])
%
% % Add intersections and comments to figure(1).
% plot(x,y,'o','MarkerEdgeColor',[0 0.6 0],'LineWidth',1)
% nuc = [2 37 39 51 68 104 106 108 110 234];     % Duplicate intersection indices.
% plot(x(nuc),y(nuc),'o','MarkerEdgeColor','k','LineWidth',1)
% text(-78,20,'Various double lines')
% text(-77,15,'Spacing')
% text(-66,15,'Crossing angle')
% title('Detection of many intersections between two curve groups')
% legend('Curve 1','Curve 2','Vertex 1','Vertex 2', ...
%                     'Intersection','Duplicate intersection','Location','SouthEast')
% ax_fig1 = axis; % Get data to match the scale of figure(2).
%
% % Extract only the regions around uncertain cases and display them.
% figure(2);     % Figure: uncertain intersection decisions only.
% nnn = size(uc_code,2);
% for n=1:nnn
%   n1 = uc_code(1,n);
%   n2 = uc_code(2,n);;
%   h1 = plot([x1(n1) x1(n1+1)],[y1(n1) y1(n1+1)],'b');
%   hold on
%   h2 = plot([x2(n2) x2(n2+1)],[y2(n2) y2(n2+1)],'r');
%   if uc_code(5,n)~=0
%     cp = uc_code(5,n);
%     h3 = plot(x(cp),y(cp),'o','MarkerEdgeColor',[0 0.6 0], ...
%                                                   'LineWidth',1);
%     text(x(cp),y(cp)+n*0.2,num2str(cp))
%   else
%     h4 = plot((x1(n1)+x1(n1+1))/2,(y1(n1)+y1(n1+1))/2,'xk', ...
%                                                   'LineWidth',1);
%     text((x1(n1)+x1(n1+1))/2,(y1(n1)+y1(n1+1))/2+n*0.2, ...
%                                       ['line(' num2str(n1) ')'])
%   end
% end
% axis equal;
% axis([ax_fig1(1) ax_fig1(2) ax_fig1(3) ax_fig1(4)]);
% title('Only regions where intersection decisions were uncertain')
% legend([h1,h2,h3,h4],'Curve 1','Curve 2','Detected intersection', ...
%                          'Uncertain curve-1 segment','Location','NorthWest')
% text(-75,-55,['Single numbers are intersection indices; ' ...
%                      'numbers after "line" are uncertain segment indices in curve 1.'])
%
% figure(1);      % Bring the main graph to the top layer.

