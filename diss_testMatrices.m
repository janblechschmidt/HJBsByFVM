% Updated 2021-11-08
% Load helper functions
clear all, clc
helperFunctions

prob.dim = 2;
prob.mode = 'MinimumArrivalTime';
% prob.mode = 'CinftySolTime';
prob.outputdir = 'output2d';
prob.trivialControl = 0;
prob.domain = 'UnitSquare';
% prob.mesh_type = 'irregular';
prob.mesh_type = 'regular';
% prob.dual_mesh = 'voronoi';
prob.dual_mesh = 'centroid';

opt.upwinding = 1;
prob.DirichletBoundary = @(onBoundary,x) onBoundary;
prob.hmax = 0.5; % Coarse initial grid
% prob.hmax = 0.25; % Coarse initial grid
% prob.hmax = 0.1;
% prob.hmax = 0.13; % Fine initial grid
% prob.hmax = 0.09; % Fine initial grid
% prob.hmax = 0.03; % Fine initial grid
% prob.grid_type = 'right';
prob.grid_type = 'left'; % positive correlation
prob.xmin = 0;
prob.xmax = 1;
prob.ymin = 0;
prob.ymax = 1;

prob = setupParameters(prob);
coef = setupCoefficients(prob);
grid = setupGrid(prob);

% The following is now part of setupGrid.m

% xlowbnd = grid.p(1,:)<=prob.xmin+1e-10;
% xupbnd = grid.p(1,:)>=prob.xmax-1e-10;
% ylowbnd = grid.p(2,:)<=prob.ymin+1e-10;
% yupbnd = grid.p(2,:)>=prob.ymax-1e-10;

% grid.D2{1}(xlowbnd | xupbnd,:) = 0;
% grid.D2{2}(xlowbnd | xupbnd | ylowbnd | yupbnd,:) = 0;
% grid.D2{3}(xlowbnd | xupbnd | ylowbnd | yupbnd,:) = 0;
% grid.D2{4}(ylowbnd | yupbnd,:) = 0;

% full(grid.D2{1}) 
% full(grid.D2{2})
% full(grid.D2{3})
% full(grid.D2{4})

% grid.D1{1}(xlowbnd | xupbnd,:) = 0;
% grid.D1{2}(ylowbnd | yupbnd,:) = 0;
% grid.D1_up{1}(xlowbnd | xupbnd,:) = 0;
% grid.D1_up{2}(ylowbnd | yupbnd,:) = 0;

% full(grid.D1{1})
% full(grid.D1{2})
