function X = frogbacksub_noise(X, n1, n2)
%FROGBACKSUB Remove FROG background.
%	FROGBACKSUB(X, [X1 X2], [Y1 Y2]) This program attempts to 
%	remove the background from a FROG trace by subtracting
%	from the trace the average of frontmost X1 and backmost 
%	X2 slices, and then the average of topmost Y1 and 
%	bottommost Y2 slices.
%
%	FROGBACKSUB(X, N1, N2) = FROGBACKSUB(X, [N1 N1], [N2 N2])
%
%	FROGBACKSUB(X, N) = FROGBACKSUB(X, [N N], [N N])
%
%	FROGBACKSUB(X) = FROGBACKSUB(X, [1 1], [1 1])

%	$Revision: 1.1 $ $Date: 2006-11-11 00:15:30 $

% error(nargchk(1,3,nargin))

switch nargin
case 1
	x1 = 1;
    x2 = 1;
    y1 = 1;
    y2 = 1;
case 2
    x1 = n1;
    x2 = n1;
    y1 = n1;
    y2 = n1;
case 3
    switch length(n1)
    case 1
        x1 = n1;
        x2 = n1;
    case 2
        x1 = n1(1);
        x2 = n1(2);
    otherwise
        error('Wrong length of input argument');
    end
    switch length(n2)
    case 1
        y1 = n2;
        y2 = n2;
    case 2
        y1 = n2(1);
        y2 = n2(2);
    otherwise
        error('Wrong length of input argument');
    end
end

[end1, end2] = size(X);
C1 = [];
C2 = [];
R1 = [];
R2 = [];

if x1 > 0
    C1 = X(:, 1 : x1);
end
if x2 > 0
    C2 = X(:, end2 - x2 + 1 : end2);
end
C = [C1, C2];
if ~isempty(C)
    Y = mean(C(:));
    X = X - Y * ones(1, end2);
end
    
if y1 > 0
    R1 = X(1 : y1, :);
end
if y2 > 0
    R2 = X(end1 - y2 + 1 : end1, :);
end
R = [R1; R2];
if ~isempty(R)
    Y = mean(R(:));
    X = X - ones(end1, 1) * Y;
end