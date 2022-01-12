function [grid] = setupZerothOrderMatrix(grid)

% Define index vectors
ix = (1:grid.N)';
jx = (1:grid.N)';

grid.D0 = sparse(ix, jx, grid.muOmega, grid.N, grid.N);

end % function [pot_mat] = setupZerothOrderMatrix(grid)
