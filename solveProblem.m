% Load helper functions
clear all%, clc
helperFunctions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solution of a Hamilton-Jacobi-Bellman equation using a Finite Volume Approach %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script implements the Finite Volume Approach for linear and particular
% non-linear PDEs of second order.
%
% This script solves one of the following problems:
%   1) elliptic non-variariatonal PDEs
%
%            [L v](x) + f(x) = 0    in Omega
%                       v(x) = g(x) on Gamma
%
%   2) parabolic non-variariatonal PDEs (backwards in time)
%
%         v_t(t,x)  + [L v](t,x) + f(t,x) = 0      in [tmin,tmax] x Omega
%                                  v(t,x) = g(t,x) on [tmin,tmax) x Gamma
%                                  v(T,x) = g(T,x) on {tmax} x Omega
%
%   3) elliptic Hamilton-Jacobi-Bellman equations
%
%         min_{u in A} { [L^u v](x) + f^u(x) } = 0    in Omega
%                                         v(x) = g(x) on Gamma
%
%   4) parabolic non-variariatonal PDEs (backwards in time)
%
%         -v_t(t,x)  - min_{u in A} { [L^u v](t,x) + f^u(t,x) } = 0      in [tmin,tmax] x Omega
%                                                        v(t,x) = g(t,x) on [tmin,tmax) x Gamma
%                                                        v(T,x) = g(T,x) on {tmax} x Omega
% 
% 
% The operator L is the infinitesimal generator of an Ito-diffusion process and given by
% [L v](x) = A(x) : D^2 v(x) + B(x)' * grad v(x) + c * v(x)
%
% Note that the discount factor c remains constant for 
opt = standardOptions();
opt.upwinding = 1;
opt.plot_solution = 1;
opt.plot_mesh = 0;
opt.plotGradient = 0;
% opt.refinementType = 'global';
opt.refinementType = 'local';
opt.number_of_refinements = 8;
opt.dual_mesh = 'voronoi';
% opt.dual_mesh = 'centroid';

% Initial mesh size
% prob.hmax = 0.02;
prob.hmax = 0.33;
% prob.hmax = 0.10;
% prob.hmax = 0.12;
% prob.hmax = 0.05;
% prob.hmax = 0.025;
% prob.hmax = 0.0125;
% prob.hmax = 0.01;

% Initial time step size (only relevant for parabolic problems)
prob.dtmax = 0.33;
% prob.dtmax = prob.hmax;
% prob.dtmax = prob.hmax^2;

prob.dual_mesh = opt.dual_mesh;

prob.mesh_type = 'regular';
% prob.mesh_type = 'irregular';

%% -----------------------------------------
%%     Problems in One Spatial Dimension
%% -----------------------------------------

%% 1) 1D: Non-variational elliptic PDES
% [prob, coef, grid] = setupNonvariationalPDE(prob, 'Poisson1d');

%% -----------------------------------------
%%     Problems in Two Spatial Dimensions
%% -----------------------------------------

%% 1) Non-variational elliptic PDES
[prob, coef, grid] = setupNonvariationalPDE(prob, 'Poisson2d');
% [prob, coef, grid] = setupNonvariationalPDE(prob, 'ConvectionDominated');
%[prob, coef, grid] = setupNonvariationalPDE(prob, 'SolInHalpha');
% [prob, coef, grid] = setupNonvariationalPDE(prob, 'CinftySol');
% [prob, coef, grid] = setupNonvariationalPDE(prob, 'NonzeroBoundary');
% [prob, coef, grid] = setupNonvariationalPDE(prob, 'LognormalDist');

%% 2) Non-variational parabolic PDES
% [prob, coef, grid] = setupParabolicNVPDE(prob, 'CinftySolTime');
[prob, coef, grid] = setupParabolicNVPDE(prob, 'WorstOfTwoAssetsPut');
 % [prob, coef, grid] = setupParabolicNVPDE(prob, 'AsianCall');

%% 3) Elliptic HJB equations
% [prob, coef, grid] = setupEllipticHJB(prob, 'MinimumArrivalTime');
% [prob, coef, grid] = setupEllipticHJB(prob, 'MinimumArrivalTimex');
% [prob, coef, grid] = setupEllipticHJB(prob, 'CinftyHJB');
% [prob, coef, grid] = setupEllipticHJB(prob, 'CinftyHJB2');

%% 4) Parabolic HJB equations

%% 5) Elliptic Variational Inequality
% [prob, coef, grid] = setupEllipticVI(prob, 'PoissonObstacle');

%% 6) Parabolic Variational Inequality
% [prob, coef, grid] = setupParabolicVI(prob, 'AmericanWorstOfTwoAssetsPut');
% [prob, coef, grid] = setupParabolicVI(prob, 'OptimalPortfolio2d');


mode = 'solve';
%mode = 'convergence';

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
return
% 
