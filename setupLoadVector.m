function rhs = setupLoadVector(grid, coef, opt, u)
if isstruct(coef)
    if nargin < 4
        if isfield(grid, 'time')
            vals = coef.f(grid.time, grid.p');
        else
            vals = coef.f(grid.p');
        end
    else
        if isfield(grid, 'time')
            vals = coef.f(grid.time, grid.p', u);
        else
            vals = coef.f(grid.p', u);
        end
    end
else
    vals = coef;
end

rhs = grid.muOmega .* vals;
