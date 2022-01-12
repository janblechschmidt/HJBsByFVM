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
opt.plot_mesh = 0;
opt.plotGradient = 0;
opt.refinementType = 'global';
opt.number_of_refinements = 8;
opt.dual_mesh = 'centroid';
prob.dual_mesh = opt.dual_mesh;

prob.mesh_type = 'irregular';
% prob.hmax = 2/3;
% prob.dtmax = 1/10;
% prob.hmax = 0.10;


% prob.hmax = 0.12;
% prob.hmax = 0.05;
% prob.hmax = 0.025;
% prob.hmax = 0.0125;
% prob.hmax = 0.01;

% Initial time step size (only relevant for parabolic problems)
% Hset = [0.1, 0.05, 0.025, 0.0125];
Hset = [0.1, 0.05, 0.025, 0.0125, 0.0125/2, 0.0125/4];
prob.dtmax = 0.001;
% prob.dtmax = prob.hmax;
% prob.dtmax = prob.hmax^2;

% Example 1
% prob.alpha = 0.0;
% prob.beta = 0.0;
% prob.sigmax = 1.0;
% prob.sigmay = 1.0;
% max(v) = 0.358438198130563

% Example 2
% prob.alpha = 0.0;
% prob.beta = 0.0;
% prob.sigmax = 0.0;
% prob.sigmay = 0.0;

% Example 3
prob.alpha = 0.2;
prob.beta = 0.2;
prob.sigmax = 0.1;
prob.sigmay = 0.1;
iter = 1;
for hmax = Hset
    if iter == 1
        prob.hmax = hmax;
        [prob, coef, grid] = setupParabolicHJB(prob, 'MATCircle');
    else
        markedTriangles = true(grid.N_tri,1);
        grid = refineIrregularGrid(prob, grid, find(markedTriangles));
    end
    
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
    iter = iter + 1;
end
