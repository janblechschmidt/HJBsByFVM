function grid = refineGrid(prob, grid, markedTriangles, markedPoints)

p = grid.p;
grid_prev = grid;

% Global refinement
if nargin < 3 | ( isempty(markedTriangles) & isempty(markedPoints) )
	switch prob.mesh_type

    case 'regular'
        % Determine previously used mesh_N (number of intervals
        % in one coordinate direction) and double its number, i.e.,
        % half the mesh size
        switch grid.dim
        case 1
            mesh_N = 2*(grid.N-1);
        case 2
            mesh_N = 2 * (sqrt(grid.N)-1);
        end
        grid = setupUniformGrid(prob, mesh_N);
        grid.c2f = setupFineToCoarseMatrix(grid, grid_prev);

    case 'irregular'
        % TODO refineIrregularGrid is not finished, it requires some other routines to be fixed.
        % In particular, the matlab function to generate a Voronoi diagram is not robust, i.e., it
        % generates duplicate entries in V, e.g., with an initial hmax of 
        % grid = refineIrregularGrid(prob, grid);
        markedTriangles=(1:grid.N_tri)';
        U = speye(grid.N,grid.N);
		[p,e,t,Unew] = refinemesh(prob.geom, grid.p, grid.e, grid.t, U,markedTriangles);
        grid.c2f = Unew;
		grid.p = p;

        % Resort edges
		grid.e = zeros(7,0);
        for i=1:size(prob.geom,2)
            ei = e(:,e(5,:)==i);
            [~,idx] = sort(ei(3,:));
            grid.e = [grid.e,ei(:,idx)];
        end

		grid.t = t;

        grid = setupIrregularGrid(prob, grid);

    end
end

grid.dim = prob.dim;
