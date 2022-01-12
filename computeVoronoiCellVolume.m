function vol = computeVoronoiCellVolume(grid)

vol = zeros(grid.N, 1);

for i=1:grid.N
	cell = grid.C{i};
	vol(i) = polyarea(grid.V(cell,1), grid.V(cell,2));
end
