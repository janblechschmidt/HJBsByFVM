clear all
helperFunctions

set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontSize',18)
% set(0,'DefaultFigureVisible','off')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Note that the discount factor c remains constant for 
opt = standardOptions();
opt.upwinding = 1;
opt.plot_solution = 1;
opt.plot_mesh = 1;
opt.plotGradient = 0;
opt.refinementType = 'global';
% opt.refinementType = 'local';
opt.number_of_refinements = 8;
% opt.dual_mesh = 'voronoi';
opt.dual_mesh = 'centroid';
prob.dual_mesh = opt.dual_mesh;
% prob.grid_type = 'right';
prob.grid_type = 'left'; % positive correlation

% prob.mesh_type = 'regular';
prob.mesh_type = 'irregular';
Hset = [0.1];
prob.dtmax = 0.05;

for hmax = Hset
    prob.hmax = hmax;
    % [prob, coef, grid] = setupEllipticHJB(prob, 'MinimumArrivalTime');
    [prob, coef, grid] = setupParabolicHJB(prob, 'OptimalPortfolio2d');
    
    mode = 'solve';
    % mode = 'convergence';
    
    switch mode
    case 'solve'
        % We plot the mesh, if opt.plot_mesh == True
        plotMesh(grid, opt);
        if strfind(prob.class, 'elliptic')
            if strfind(prob.class, 'HJB')
                sol = solveEllipticHJB(prob, grid, coef, opt);
            else
                if strfind(prob.class, 'VI')
                    sol = solveEllipticVI(prob, grid, coef, opt);
                else
                    sol = solveNVPDE(prob, grid, coef, opt);
                end
            end
        else
            sol = solveParabolicProblem(prob, grid, coef, opt);
        end
    
        % We plot the results, if opt.plot_solution == True 
        plotSolution(prob, grid, coef, opt, sol);
    
    case 'convergence'
        [sol, conv, grid] = determineConvergenceNVPDE(prob, grid, coef, opt);
        conv
    end % switch mode
    % interpolate([100,100],sol.v,grid)
    % 
    % if isfield(coef, 'solution')
    %     fprintf('Analytical solution: \n')
    %     if strfind(prob.class, 'elliptic')
    %         x = [0.5,0.5];
    %         coef.solution(x)
    %     else
    %         x = [40,40];
    %         t = prob.Tmin;
    %         coef.solution(t,x)
    %     end
    %     fprintf('Numerical solution at x: '), x
    %     fprintf('\nv(x) = \n'), interpolate(x, sol.v, grid)
    % end
    % x = [40,40];
    % t = prob.Tmin;
    % fprintf('Numerical solution at x: '), x
    % fprintf('\nv(x) = \n'), interpolate(x, sol.v, grid)
    
    % 
    % n = sqrt(size(grid.p,2))-1
    % x=linspace(prob.xmin, prob.xmax, n + 1);
    % y=linspace(prob.ymin, prob.ymax, n + 1);
    % [X,Y] = meshgrid(x,y);
    % vmat = reshape(v,n+1,n+1);
    % sfigure(10), clf
    % contour(X,Y,vmat,50)
    fprintf('Grid.N: %d\n', grid.N);
    fprintf('Check: %8.6f\n', abs(max(sol.v) - 0.78648289));
    fprintf('Max: %8.6f\n', max(sol.v));
    fprintf('Diff to 1: %8.6f\n', 1-max(sol.v));
    keyboard
end
