clear all
helperFunctions

% Note that the discount factor c remains constant for 
opt = standardOptions();
opt.upwinding = 1;
opt.plot_solution = 0;
opt.plot_mesh = 0;
opt.plotGradient = 0;
opt.refinementType = 'global';
% opt.refinementType = 'local';
opt.number_of_refinements = 8;


Hset = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.0125/2];
% Hset = [0.2];
% prob.dtmax = 0.001;
% Hset = [0.2, 0.1, 0.05, 0.025, 0.0125];
% prob.dtmax = 0.2;
% Test 1
prob.alpha = 0.0;
prob.beta = 0.0;
prob.sigmax = 0.0;
prob.sigmay = 0.0;
p.corrxy = 0.0;
vqoiexact = 1.0;
fnamebase = 'MATRhoZero';

% Test 2
% run 1
% prob.alpha = 0.5;
% prob.beta = 0.5;
% prob.sigmax = 0.5;
% prob.sigmay = 0.5;
% p.corrxy = -0.9;
% vqoiexact = 0.957254798798650;
% fnamebase = 'MATRhoNeg';

% Test 3
% run 3
% prob.alpha = 0.5;
% prob.beta = 0.5;
% prob.sigmax = 0.5;
% prob.sigmay = 0.5;
% prob.corrxy = 0.9;
% vqoiexact = 0.957254798798650;
% fnamebase = 'MATRhoPos';
% load('../hjb_finite_difference/fine_sol_MAT_Pos.mat')
% 1 run
% 2 run
% 3 run
% 4 run
% 5 open


% Test 4: non-constant correlation
% run 2
prob.alpha = 0.5;
prob.beta = 0.5;
prob.sigmax = 0.5;
prob.sigmay = 0.5;
prob.corrxy = 0.9;
fnamebase = 'MATRhoSwitch';
load('../hjb_finite_difference/fine_sol_MAT_Switch.mat')
% 1
% 2 
% 3 open
% 4
% 5 open

fd_n = sqrt(length(fd_v));
fd_vmat = reshape(fd_v,fd_n,fd_n);
fd_X = reshape(fd_XY(:,1),fd_n,fd_n);
fd_Y = reshape(fd_XY(:,2),fd_n,fd_n);
vqoiexact = interp2(fd_X,fd_Y,fd_vmat,[0],[0]);

for pmode = [5]
% pmode = 2;
switch pmode
    case 1
    prob.mesh_type = 'regular';
    prob.grid_type = 'left'; % positive correlation
    opt.dual_mesh = 'centroid';
    fname = [fnamebase, '_UniLeftCent'];
    case 2
    prob.mesh_type = 'regular';
    prob.grid_type = 'left'; % positive correlation
    opt.dual_mesh = 'voronoi';
    fname = [fnamebase, '_UniLeftVoro'];
    case 3
    prob.mesh_type = 'regular';
    prob.grid_type = 'right'; % positive correlation
    opt.dual_mesh = 'centroid';
    fname = [fnamebase, '_UniRightCent'];
    case 4
    prob.mesh_type = 'regular';
    prob.grid_type = 'right'; % positive correlation
    opt.dual_mesh = 'voronoi';
    fname = [fnamebase, '_UniRightVoro'];
    case 5
    prob.mesh_type = 'irregular';
    opt.dual_mesh = 'centroid';
    fname = [fnamebase, '_IrrCent'];
end
prob.dual_mesh = opt.dual_mesh;

% prob.hmax = 0.12;
% prob.hmax = 0.05;
% prob.hmax = 0.025;
% prob.hmax = 0.0125;
% prob.hmax = 0.01;

% Initial time step size (only relevant for parabolic problems)
% Hset = [2/200];
% prob.dtmax = 0.01;
% prob.dtmax = prob.hmax;
% prob.dtmax = prob.hmax^2;

% Example 1
% Hset = [2/200];


% irregular voronoi: 0.786704
% irregular centroid: 0.761544
% regular voronoi: 0.786483
% regular centroid: 0.711387
iter = 1;
for hmax = Hset
    fprintf('Next hmax\n\n');
    if strcmp(prob.mesh_type,'regular') | iter == 1
        prob.hmax = hmax;
        prob.dtmax = hmax/10;
        [prob, coef, grid] = setupParabolicHJB(prob, 'ParabolicMinimumArrivalTime');
    else
        markedTriangles = true(grid.N_tri,1);
        grid = refineIrregularGrid(prob, grid, find(markedTriangles));
        prob.dtmax = hmax/10;
    end
    
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
    % fprintf('Check: %8.6f\n', abs(max(sol.v) - 0.78648289));
    fprintf('Max: %8.6f\n', max(sol.v));
    fprintf('Diff to 1: %8.6f\n', 1-max(sol.v));
    vglob = zeros(fd_n^2,1);
    for j=1:fd_n
        fprintf('j/fd_n: %d/%d\n', j, fd_n);
        jidx = ((j-1)*fd_n+1):(j*fd_n);
        vglob(jidx) = interpolate(fd_XY(jidx,:), sol.v, grid);
    end
    vqoi = interpolate([0.0,0.0],sol.v,grid);
    out(iter,:) = [iter, grid.N, sol.nt, max(sol.v), max(abs(vglob-fd_v)), vqoiexact-max(sol.v)];
    out(:,5)
    % keyboard
    iter = iter + 1;
end

fheader={'Iter', 'Ndofs', 'Nt', 'max_v', 'e_comp','e_qoi'};
fname = sprintf('%supwinding_%d_%d', fname, opt.upwinding, grid.N);
fname = strcat('./tables/', fname, '.csv');
fid = fopen(fname, 'w');
fprintf(fid, '%s\n', strjoin(fheader,','));
dlmwrite(fname, out, '-append','precision',10);
end
