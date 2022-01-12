function C = getZerothOrderMatrix(prob, grid, coef, u)

if nargin < 4
    u = [];
end
if isstruct(coef)
    vals = evaluateCoefficient(prob, grid, coef.potential, u);
else
    vals = coef;
end
C = spdiags(vals, [0], grid.N, grid.N) * grid.D0;

end % function C = getZerothOrderMatrix(prob, grid, coef, u)
