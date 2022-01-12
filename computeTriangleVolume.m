function [V, orientation] = computeTriangleVolume(grid)
% Compute volumes of triangles in primal grid

p = grid.p;
t = grid.t';
x1 = p(1, t(:,1))';
x2 = p(1, t(:,2))';
x3 = p(1, t(:,3))';
y1 = p(2, t(:,1))';
y2 = p(2, t(:,2))';
y3 = p(2, t(:,3))';

V = 0.5 * ((x3-x1) .* (y2-y1) - (x2-x1) .* (y3-y1));
if nargout == 2
	orientation = zeros(size(V));
	orientation(V > 0) = -1;
	orientation(V < 0) = 1;
end
V = abs(V);

end % function V = computeTriangleVolume(grid)
