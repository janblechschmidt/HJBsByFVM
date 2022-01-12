function [sol] = solveParabolicProblem(prob, grid, coef, opt, sol)

% Initialize solution vector with correct boundary and final time conditions

t = prob.Tmax;

sol.v = coef.g(prob.Tmax, grid.p');
sol.nt = ceil((prob.Tmax-prob.Tmin)/prob.dtmax)+1;
Tspace = prob.Tmax - linspace(0, prob.Tmax-prob.Tmin, sol.nt);

for i = 1:length(Tspace)-1
    t = Tspace(i+1);
    fprintf('Time step %3d (t = %4.2f)\n', i, t)

    dt = Tspace(i) - Tspace(i+1);
    grid.dt = dt;
    grid.time = t;
    if strfind(prob.class, 'HJB')
        sol = solveEllipticHJB(prob, grid, coef, opt, sol);
    else
        if strfind(prob.class, 'VI')
            sol = solveEllipticVI(prob, grid, coef, opt, sol);
        else
            sol = solveNVPDE(prob, grid, coef, opt, sol);
        end
    end

end % for i = 1:length(Tspace)-1
end % function [v] = solveParabolicProblem(prob, grid, coef, opt)
