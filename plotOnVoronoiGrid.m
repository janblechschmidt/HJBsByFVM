function ax = plotOnVoronoiGrid(grid, f, ax)
% Plot function f defined on Voronoi cells
hold on
for i = 1:grid.N
	cell = grid.C{i};
	X = grid.V([cell,cell(1)],1);
	Y = grid.V([cell,cell(1)],2);
	Z = f(i) * ones(size(X));
	C = Z;
	fill3(X,Y,Z,C,'EdgeAlpha',0.1,'LineWidth',0.01);
end

% function ax = plotOnVoronoiGrid(grid, f, ax)
