function A = getSecondOrderMatrix(prob, grid, coef, u)

if nargin < 4
    u = [];
end

if isstruct(coef)
    vals = evaluateCoefficient(prob, grid, coef.diffusion, u);
else
    vals{1} = coef(:,1);
    vals{2} = coef(:,2);
    vals{3} = coef(:,3);
    vals{4} = coef(:,4);
end

A = sparse(grid.N,grid.N);
for m = 1:grid.dim
    for n = 1:grid.dim
        % Compute derivative v_{x_m,x_n}
        idx = grid.dim*(m-1)+n;
        % TODO Check, whether vals has to be indexed according to triangle enumeration or not,
        % previously vij = vij(i) with vij = vals{idx};
        Aij = spdiags(grid.muOmega.*vals{idx}, [0], grid.N, grid.N) * grid.D2{idx};
        A = A + Aij;
    end % for j = 1:grid.dim
end % for i = 1:grid.dim

end % function A = getSecondOrderMatrix(prob, grid, coef, u)
