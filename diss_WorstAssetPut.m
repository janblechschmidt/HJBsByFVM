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
opt.plot_solution = 0;
opt.plot_mesh = 0;
opt.plotGradient = 0;
opt.refinementType = 'global';
% opt.refinementType = 'local';
opt.number_of_refinements = 8;
% opt.dual_mesh = 'voronoi';
opt.dual_mesh = 'centroid';
prob.dual_mesh = opt.dual_mesh;
% prob.grid_type = 'right'; % negative correlation
prob.grid_type = 'left'; % positive correlation

% prob.mesh_type = 'regular';
% prob.mesh_type = 'irregular';

mode = 2;
switch mode
    case 1
    % First setup (equivalent to fd scheme)
    opt.dual_mesh = 'voronoi';
    prob.dual_mesh = opt.dual_mesh;
    prob.grid_type = 'left'; % positive correlation
    prob.mesh_type = 'regular';
    % prob.mesh_type = 'irregular';
    Hset = [20,10,5,2.5,1.25,1.25/2, 1.25/4];
    prob.dtmax = 0.0005;
    fnamebase = 'WorstAssetUniformVoronoi_'
    custom_refine = false;

    case 2
    % First setup (equivalent to fd scheme)
    opt.dual_mesh = 'centroid';
    prob.dual_mesh = opt.dual_mesh;
    prob.grid_type = 'left'; % positive correlation
    prob.mesh_type = 'regular';
    % prob.mesh_type = 'irregular';
    % Hset = [20,10,5,2.5,1.25,1.25/2];
    Hset = [20,10,5,2.5,1.25,1.25/2, 1.25/4];
    prob.dtmax = 0.0005;
    fnamebase = 'WorstAssetUniformCentroid_'
    custom_refine = false;

    case 3
    % Second setup (locally refined)
    opt.dual_mesh = 'centroid';
    prob.dual_mesh = opt.dual_mesh;
    % prob.grid_type = 'left'; % positive correlation
    % prob.mesh_type = 'regular';
    prob.mesh_type = 'irregular';
    Hset = [40,20,10,5,2.5,1.25];
    prob.dtmax = 0.0005;
    fnamebase = 'WorstAssetRefinedCentroid_'
    custom_refine = true;
end

iter = 1;

for hmax = Hset
    if strcmp(prob.mesh_type,'regular') | iter == 1
        prob.hmax = hmax;
        [prob, coef, grid] = setupParabolicNVPDE(prob, 'WorstOfTwoAssetsPut');
        % Custom refinement
        if custom_refine
            % Triangle centers
            tc = 1/3*(grid.p(:,grid.t(1,:))+ grid.p(:,grid.t(2,:))+grid.p(:,grid.t(3,:)));
            markedTriangles = false(grid.N_tri,1);
            markedTriangles(tc(1,:) > 20 & tc(1,:) < 60 & tc(2,:) > 20) = true;
            markedTriangles(tc(2,:) > 20 & tc(2,:) < 60 & tc(1,:) > 20) = true;
            markedTriangles(abs(tc(2,:)-tc(1,:)) < 30 & tc(2,:) < 60 & tc(1,:) < 60) = true;
            grid = refineIrregularGrid(prob, grid, find(markedTriangles));
            tc = 1/3*(grid.p(:,grid.t(1,:))+ grid.p(:,grid.t(2,:))+grid.p(:,grid.t(3,:)));
            markedTriangles = false(grid.N_tri,1);
            markedTriangles(tc(1,:) > 30 & tc(1,:) < 50 & tc(2,:) > 30) = true;
            markedTriangles(tc(2,:) > 30 & tc(2,:) < 50 & tc(1,:) > 30) = true;
            markedTriangles(abs(tc(2,:)-tc(1,:)) < 15 & tc(2,:) < 50 & tc(1,:) < 50) = true;
            grid = refineIrregularGrid(prob, grid, find(markedTriangles));
            plotMesh(grid, opt);
        end
    else
        markedTriangles = true(grid.N_tri,1);
        grid = refineIrregularGrid(prob, grid, find(markedTriangles));
    end

    fprintf('start solve')
    % return
    
    mode = 'solve';
    % mode = 'convergence';
    
    switch mode
    case 'solve'
        % We plot the mesh, if opt.plot_mesh == True
        % plotMesh(grid, opt);
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
        vqoi = interpolate([40,40],sol.v,grid);
        fprintf('Value at qoi: %8.4f\n',vqoi);
        % We plot the results, if opt.plot_solution == True 
        plotSolution(prob, grid, coef, opt, sol);
        NN = 160;
        [X,Y] = meshgrid(linspace(30,50,NN+1),linspace(30,50,NN+1));
        XY = [reshape(X,[],1),reshape(Y,[],1)];
        Vint = interpolate(XY, sol.v, grid)';
        Vex = coef.solution(prob.Tmin,XY);
        figure(5); clf;
        surf(X,Y,reshape(abs(Vex-Vint),NN+1,NN+1),'EdgeAlpha',0.1);
        xlabel('$x$')
        ylabel('$y$')
        colormap viridis
        view(2)
        colorbar('eastoutside')
        axis equal
        axis tight
        caxis([0.01,0.04])
        fname = sprintf('./tex/%supwinding_%d_%d.png', fnamebase,opt.upwinding, grid.N);
        print(gcf,fname,'-dpng','-r300'); 
        vexact = coef.solution(prob.Tmin,grid.p');
        out(iter,:) = [iter, grid.N, sol.nt, max(abs(vexact-sol.v)),max(abs(Vex-Vint)), max(abs(vqoi-coef.solution(prob.Tmin,[40,40])))];
        % surf(X,Y,reshape(Vint,NN+1,NN+1));

    
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

    if isfield(coef,'solution')
        if isfield(grid, 'time')
            vexact = coef.solution(grid.time,grid.p');
        else
            if isfield(prob,'Tmax')
                vexact = coef.solution(prob.Tmin, grid.p');
            else
                vexact = coef.solution(grid.p');
            end
        end
        fprintf('Maximum error: %6.4e\n', max(abs(vexact - sol.v)));
    end

    iter = iter + 1;
    % keyboard
end
fheader={'Iter', 'Ndofs', 'Nt', 'e_global', 'e_local', 'e_qoi'};
fname = sprintf('%supwinding_%d_%d', fnamebase, opt.upwinding, grid.N);
fname = strcat('./tables/', fname, '.csv');
fid = fopen(fname, 'w');
fprintf(fid, '%s\n', strjoin(fheader,','));
dlmwrite(fname, out, '-append','precision',10);
