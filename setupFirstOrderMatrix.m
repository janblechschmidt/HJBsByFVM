function [grid] = setupFirstOrderMatrix(grid)

Minv = spdiags(1./grid.muOmega, [0], grid.N, grid.N);
% Minv = speye(grid.N);
mij = grid.muGamma;

for i = 1:grid.dim
    K = grid.Nuij{i};

    % Define central difference scheme, i.e.,
    % v_{i,j} = 0.5*(v_i + v_j)
    B_offdiag = 0.5*K;
    B_diag = 0.5*spdiags(sum(K,2), [0], grid.N, grid.N);
    Bi = B_diag + B_offdiag;
    grid.D1{i} = Minv*Bi;
    % Define matrix for first order upwind scheme, i.e.,
    % during runtime, the values are assigned either to
    % v_i or v_j
    grid.D1_up{i} = Minv*K;
    % keyboard
    tmpforw = K.*(K>0);
    grid.D1forw{i} = Minv*(tmpforw - spdiags(sum(tmpforw,2), [0], grid.N, grid.N));
    tmpback = K.*(K<0);
    grid.D1back{i} = Minv*(tmpback - spdiags(sum(tmpback,2), [0], grid.N, grid.N));
    
    % Matrix containing domain boundary information
    tb = grid.n_boundary;
    % keyboard
    % First: Use central weighting on edges adjacent to inner vertices
    % grid.D1forw{i}(grid.bnd,:) = grid.D1{i}(grid.bnd,:);
    % grid.D1forw{i}  = grid.D1forw{i} +  Minv*sparse(tb(3,:),tb(3,:),tb(i,:).*tb(4,:),grid.N,grid.N);
    % grid.D1back{i}(grid.bnd,:) = grid.D1{i}(grid.bnd,:);
    % grid.D1back{i} = grid.D1back{i} + Minv*sparse(tb(3,:),tb(3,:),tb(i,:).*tb(4,:),grid.N,grid.N);
    % Previous
    % grid.D1forw{i}(grid.lowbnd{i} | grid.upbnd{i},:) = grid.D1{i}(grid.lowbnd{i} | grid.upbnd{i},:);
    % grid.D1back{i}(grid.lowbnd{i} | grid.upbnd{i},:) = grid.D1{i}(grid.lowbnd{i} | grid.upbnd{i},:);
    grid.D1forw{i}(grid.DomainBoundary,:) = grid.D1{i}(grid.DomainBoundary,:);
    grid.D1back{i}(grid.DomainBoundary,:) = grid.D1{i}(grid.DomainBoundary,:);
    grid.Dout{i} = Minv*sparse(tb(3,:),tb(3,:),tb(i,:).*tb(4,:),grid.N,grid.N);
    % keyboard
    grid.D1forw{i} = grid.D1forw{i} + grid.Dout{i};
    grid.D1back{i} = grid.D1back{i} + grid.Dout{i};
    % keyboard
    % grid.D1forw{i} = Minv*(K.*(K>0) + spdiags(sum(K.*(K<0),2), [0], grid.N, grid.N))
    % grid.D1back{i} = Minv*(K.*(K<0) + spdiags(sum(K.*(K>0),2), [0], grid.N, grid.N))
    % grid.D1{i}(grid.bnd,:) = grid.D1{i}(grid.bnd,:) / 2;
    % grid.D1_up{i}(grid.bnd,:) = grid.D1_up{i}(grid.bnd,:) / 2;


end
% Ensure stencils are only applied where they are usefull
% keyboard
if isfield(grid, 'lowbnd')
    grid.D1{1}(grid.lowbnd{1} | grid.upbnd{1},:) = 0;
    grid.D1{2}(grid.lowbnd{2} | grid.upbnd{2},:) = 0;
    
    grid.D1_up{1}(grid.lowbnd{1} | grid.upbnd{1},:) = 0;
    grid.D1_up{2}(grid.lowbnd{2} | grid.upbnd{2},:) = 0;
end

% grid.D1forw{1}(grid.xupbnd,:) = 0;
% grid.D1forw{2}(grid.yupbnd,:) = 0;

% grid.D1back{1}(grid.xlowbnd,:) = 0;
% grid.D1back{2}(grid.ylowbnd,:) = 0;

% grid.D1forw{1}(grid.xlowbnd | grid.xupbnd,:) = grid.D1forw{1}(grid.xlowbnd | grid.xupbnd,:) / 2;
% grid.D1back{1}(grid.xlowbnd | grid.xupbnd,:) = grid.D1back{1}(grid.xlowbnd | grid.xupbnd,:) / 2;
% grid.D1forw{2}(grid.ylowbnd | grid.yupbnd,:) = grid.D1forw{2}(grid.ylowbnd | grid.yupbnd,:) / 2;
% grid.D1back{2}(grid.ylowbnd | grid.yupbnd,:) = grid.D1back{2}(grid.ylowbnd | grid.yupbnd,:) / 2;

% grid.D1forw{1}(grid.xupbnd,:) = grid.D1back{1}(grid.xupbnd,:);
% grid.D1forw{2}(grid.yupbnd,:) = grid.D1back{2}(grid.yupbnd,:);

% grid.D1back{1}(grid.xlowbnd,:) = grid.D1forw{1}(grid.xlowbnd,:);
% grid.D1back{2}(grid.ylowbnd,:) = grid.D1forw{2}(grid.ylowbnd,:);
end % function [grid] = setupFirstOrderMatrix(grid)
