function grid = setupSecondOrderMatrix(grid)
switch grid.dim
case 1
    fprintf('inside setupSecondOrderMatrix\n')
    e = ones(grid.N,1);
    scale = 0.5*([grid.muT;0] + [0;grid.muT]);
    % tmp = spdiags(1./grid.muT.*[[e(2:end);111], -2*e, [111;e(2:end)]], [-1,0,1], grid.N, grid.N)
    s = grid.muT;
    ix = [2:grid.N-1, 2:grid.N-1, 2:grid.N-1, 2:grid.N-1];
    iy = [2:grid.N-1, 3:grid.N,   2:grid.N-1, 1:grid.N-2];
    iz = [-1./grid.muT(2:end)', 1./grid.muT(2:end)', -1./grid.muT(1:end-1)', 1./grid.muT(1:end-1)'];
    
    grid.D2{1} = spdiags(1./grid.muOmega,[0],grid.N,grid.N)*sparse(ix, iy, iz, grid.N, grid.N);

case 2
    t = grid.t';
    
    T = [t(:,[1,2,3]); t(:,[2,3,1]); t(:,[3,1,2])];
    
    V = repmat(grid.muT,3,1);
    O = repmat(grid.orientationT, 3, 1);
    
    i = T(:,1);
    k = T(:,2);
    j = T(:,3);
    
    % The following step can be simplified, work is done threefold
    x_kj = grid.p(:,j)'-grid.p(:,k)';
    x_ji = grid.p(:,i)'-grid.p(:,j)';
    x_ik = grid.p(:,k)'-grid.p(:,i)';
    % x_kj + x_ji + x_ik % should be zero 
    
    % Ensure, that ni is outer normal vector and not inner normal vector
    li = sqrt(sum(x_kj.^2,2));
    ni = ([-x_kj(:,2), x_kj(:,1)] ./ li) .*O;
    
    lk = sqrt(sum(x_ji.^2,2));
    nk = ([-x_ji(:,2), x_ji(:,1)] ./ lk) .*O;
    
    lj = sqrt(sum(x_ik.^2,2));
    nj = ([-x_ik(:,2), x_ik(:,1)] ./ lj) .*O;
    % ni .* li + nj .* lj + nk .* lk % should be zero 
    
    % Initialize matrix for the second order term
    Ncon = grid.N_tri * 3;
    for m = 1:grid.dim
        for n = 1:grid.dim
            % Compute derivative v_{x_m,x_n}
            idx = grid.dim*(m-1)+n;
            valij = - 1./(4*V) .* li .* lj .* ( ni(:,m) .* ones(grid.N_tri * 3,1) .* nj(:,n) );
            valik = - 1./(4*V) .* li .* lk .* ( ni(:,m) .* ones(grid.N_tri * 3,1) .* nk(:,n) );
            tmp = sparse([i;i], [j;k], [valij; valik], grid.N, grid.N);
            tmp = tmp - spdiags(sum(tmp,2), [0], grid.N, grid.N);
            grid.D2{idx} = spdiags(1./grid.muOmega,[0],grid.N,grid.N)*tmp;
        end % for j = 1:grid.dim
    end % for i = 1:grid.dim
end % switch grid.dim

% Set rows to zero along boundaries
% grid.D2{1}(grid.lowbnd{1} | grid.upbnd{1},:) = 0;
% grid.D2{2}(grid.lowbnd{1} | grid.upbnd{1} | grid.lowbnd{2} | grid.upbnd{2},:) = 0;
% grid.D2{3}(grid.lowbnd{1} | grid.upbnd{1} | grid.lowbnd{2} | grid.upbnd{2},:) = 0;
% grid.D2{4}(grid.lowbnd{2} | grid.upbnd{2},:) = 0;

end % function
