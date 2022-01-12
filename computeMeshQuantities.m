function grid = computeMeshQuantities(grid)

% Compute volumes of triangles
[grid.muT, grid.orientationT] = computeTriangleVolume(grid);

grid = computeDualMeshQuantities(grid);
% 
% % Compute volumes of Voronoi cells
% grid.muOmega = computeVoronoiCellVolume(grid);
% 
% % Compute lengths of edges between Voronoi cells
% grid.muGamma = computeVoronoiEdgeLength(grid);
% 
% % Connectivity for primal grid
% for i=1:grid.N
% 	Ni{i} = find(grid.muGamma(i,:)>0);
% end
% grid.Ni = Ni;
% 
% % Determine some further auxiliary matrices
% eN = ones(grid.N,1);
% % D_i,j holds the distance between points i and j
% D = sqrt((grid.p(1,:)-grid.p(1,:)').^2+(grid.p(2,:)-grid.p(2,:)').^2);
% 
% % Matrices holding the x and y coordinates of midpoints along
% % triangles
% % TODO: Check, whether these matrices are necessary any more
% % grid.Mx = .5 * (grid.p(1,:)' * eN' + eN * grid.p(1,:));
% % grid.My = .5 * (grid.p(2,:)' * eN' + eN * grid.p(2,:));
% 
% % Matrices holding the x and y coordinates of normalized
% % vectors pointing from x_i to x_j
% grid.Nij{1} = (grid.p(1,:)-grid.p(1,:)') ./ D;
% grid.Nij{2} = (grid.p(2,:)-grid.p(2,:)') ./ D;
% 
% % Replace the assignment of NaN on the diagonals by zero
% grid.Nij{1}(1:grid.N+1:end) = 0;
% grid.Nij{2}(1:grid.N+1:end) = 0;
% 
end
