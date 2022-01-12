function vals = evaluateCoefficient(prob, grid, fun, u)
if iscell(fun)
    if strfind(prob.class, 'elliptic')
        if strfind(prob.class, 'HJB')
            vals = cellfun(@(f) f(grid.p', u), fun, 'UniformOutput', false);
        else
            vals = cellfun(@(f) f(grid.p'), fun, 'UniformOutput', false);
        end
    else
        if strfind(prob.class, 'HJB')
            vals = cellfun(@(f) f(grid.time, grid.p', u), fun, 'UniformOutput', false);
        else
            vals = cellfun(@(f) f(grid.time, grid.p'), fun, 'UniformOutput', false);
        end
    end
else
    if strfind(prob.class, 'elliptic')
        if strfind(prob.class, 'HJB')
            vals = fun(grid.p', u);
        else
            vals = fun(grid.p');
        end
    else
        if strfind(prob.class, 'HJB')
            vals = fun(grid.time, grid.p', u);
        else
            vals = fun(grid.time, grid.p');
        end
    end
end
