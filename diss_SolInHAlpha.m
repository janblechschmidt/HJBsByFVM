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
opt.upwinding = 0;
opt.plot_solution = 0;
opt.plotGradient = 0;
opt.refinementType = 'global';
% opt.refinementType = 'local';
opt.number_of_refinements = 8;


opt.plot_mesh = 0;
mode = 3;
switch mode
    case 1
    prob.mesh_type = 'regular';
    opt.dual_mesh = 'voronoi';
    fnamebase = 'SolInHAlpha_2UniVoro'
    custom_refine=false;
    % Hset = [0.2,0.1,0.05,0.025,0.0125,0.0125/2,0.0125/4];
    Hset = [1/6.25,1/12.5,1/25,1/50,1/100,1/200,1/400];
    prob.grid_type = 'left'; % positive correlation
    case 2
    prob.mesh_type = 'regular';
    opt.dual_mesh = 'centroid';
    fnamebase = 'SolInHAlpha_2UniCent'
    custom_refine=false;
    % Hset = [0.2,0.1,0.05,0.025,0.0125,0.0125/2,0.0125/4];
    Hset = [1/6.25,1/12.5,1/25,1/50,1/100,1/200,1/400];
    % Hset = [0.2,0.1,0.05,0.025,0.0125,0.0125/2];
    prob.grid_type = 'left'; % positive correlation
    case 3
    prob.mesh_type = 'irregular';
    opt.dual_mesh = 'centroid';
    fnamebase = 'SolInHAlpha_2RefCent'
    custom_refine=true;
    Hset = [0.22,0.1,0.05,0.025,11,12,13];
end
prob.dual_mesh = opt.dual_mesh;
    
iter = 1;

for hmax = Hset
    fprintf('Current hmax %8.4f\n', hmax)
    if strcmp(prob.mesh_type,'regular') | iter == 1
        prob.hmax = hmax;
        [prob, coef, grid] = setupNonvariationalPDE(prob, 'SolInHalpha');
        % Custom refinement
        if custom_refine
            % Triangle centers
            % tc = 1/3*(grid.p(:,grid.t(1,:))+ grid.p(:,grid.t(2,:))+grid.p(:,grid.t(3,:)));
            % rc = sqrt(sum(tc.^2,1));
            % markedTriangles = false(grid.N_tri,1);
            % markedTriangles(rc < .7) = true;
            % grid = refineIrregularGrid(prob, grid, find(markedTriangles));

            tc = 1/3*(grid.p(:,grid.t(1,:))+ grid.p(:,grid.t(2,:))+grid.p(:,grid.t(3,:)));
            rc = sqrt(sum(tc.^2,1));
            markedTriangles = false(grid.N_tri,1);
            markedTriangles(rc < .75) = true;
            markedTriangles(sqrt(sum((tc-[1;1]).^2,1)) <0.45) = true;
            markedTriangles(sqrt(sum((tc-[0;1]).^2,1)) <0.45) = true;
            markedTriangles(sqrt(sum((tc-[1;0]).^2,1)) <0.45) = true;
            grid = refineIrregularGrid(prob, grid, find(markedTriangles));

            tc = 1/3*(grid.p(:,grid.t(1,:))+ grid.p(:,grid.t(2,:))+grid.p(:,grid.t(3,:)));
            rc = sqrt(sum(tc.^2,1));
            markedTriangles = false(grid.N_tri,1);
            markedTriangles(rc < .4) = true;
            grid = refineIrregularGrid(prob, grid, find(markedTriangles));

            tc = 1/3*(grid.p(:,grid.t(1,:))+ grid.p(:,grid.t(2,:))+grid.p(:,grid.t(3,:)));
            rc = sqrt(sum(tc.^2,1));
            markedTriangles = false(grid.N_tri,1);
            markedTriangles(rc < .2) = true;
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
        vqoi = interpolate([0.5,0.5],sol.v,grid);
        % fprintf('Value at qoi: %8.4f\n',vqoi);
        % We plot the results, if opt.plot_solution == True 
        plotSolution(prob, grid, coef, opt, sol);
        NN = 160;
        [X,Y] = meshgrid(linspace(0,0.2,NN+1),linspace(0,0.2,NN+1));
        XY = [reshape(X,[],1),reshape(Y,[],1)];
        Vint = interpolate(XY, sol.v, grid)';
        Vex = coef.solution(XY);
        figure(5); clf;
        surf(X,Y,reshape(abs(Vex-Vint),NN+1,NN+1),'EdgeAlpha',0.1);
        xlabel('$x$')
        ylabel('$y$')
        colormap viridis
        view(2)
        colorbar('eastoutside')
        axis equal
        axis tight
        % caxis([0,0.1])
        fname = sprintf('./tex/%supwinding_%d_%d.png', fnamebase,opt.upwinding, grid.N);
        print(gcf,fname,'-dpng','-r300'); 
        vexact = coef.solution(grid.p');
        % surf(X,Y,reshape(Vint,NN+1,NN+1));

        out(iter,:) = [iter, grid.N, max(abs(vexact-sol.v)),max(abs(Vex-Vint)), max(abs(vqoi-coef.solution([0.5,0.5])))];
    
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
fheader={'Iter', 'Ndofs', 'e_global', 'e_local', 'e_qoi'};
fname = sprintf('%supwinding_%d_%d', fnamebase, opt.upwinding, grid.N);
fname = strcat('./tables/', fname, '.csv');
fid = fopen(fname, 'w');
fprintf(fid, '%s\n', strjoin(fheader,','));
dlmwrite(fname, out, '-append','precision',10);
