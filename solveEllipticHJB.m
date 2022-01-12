function [sol] = solveEllipticHJB(prob, grid, coef, opt, sol)

% Check whether problem is parabolic or elliptic, this changes
% the behaviour of the code slightly.

if strfind(prob.class, 'parabolic')
    v = sol.v;
    v_np1 = sol.v;
    fac_dt = 1/grid.dt;
else
    v = coef.g(grid.p');
    v_np1 = zeros(size(v));
    fac_dt = 0;
end

% Initialize Howard iteration
Howard_done = 0;
Howard_max_iter = 50;
Howard_eps = 1e-8;
Howard_iter = 0;

while ~Howard_done

    u = coef.control_law(prob, grid, coef, v, opt.upwinding);
    % TODO This has to be updated since do not want to rely on the control law.
    % gradv = cellfun(@(C) C*v,grid.D1,'UniformOutput',false');
    % hessv = cellfun(@(C) C*v,grid.D2,'UniformOutput',false');
    % if isfield(grid,'time')
    %     u = coef.control_law(grid.time, grid.p', v, gradv, hessv);
    % else
    %     u = coef.control_law(grid.p', v, gradv, hessv);
    % end
    % Plot current control
    % sfigure(6); clf;
    % for i = 1:prob.control_dim
    %     subplot(prob.control_dim, 1, i);
    %     plotOnVoronoiGrid(grid,u(:,i));
    %     title('u')
    % end

    % Plot gradient
    if opt.plotGradient
        sfigure(7); clf;
        for i = 1:prob.dim
            subplot(prob.dim, 1, i);
            plotOnTriangulation(grid.p, grid.t' ,gradv{i});
            title('gradv')
        end
        pause
    end % if opt.plotGradient

    % Plot Hessian matrix
    if opt.plotHessian
        sfigure(8); clf;
        for i = 1:prob.dim
            for j = 1:prob.dim
                subplot(prob.dim, prob.dim, 2*(i-1)+j);
                % tmp = spdiags(1./grid.muOmega,[0],grid.N,grid.N)*hessv{2*(i-1)+j};
                plotOnTriangulation(grid.p, grid.t' ,hessv{2*(i-1)+j});
                title(sprintf('D^2_{x%ix%i} v', i, j));
                % keyboard
            end
        end
        pause
    end % if opt.plotHessian
    
    v_prev = v;

    % Assemble coeffcient matrices
    A = getSecondOrderMatrix(prob, grid, coef, u);
    B = getFirstOrderMatrix(prob, grid, coef, opt.upwinding, u);
    C = getZerothOrderMatrix(prob, grid, coef, u);
    T = grid.D0;
    
    % Assemble right-hand side
    rhs = - setupLoadVector(grid, coef, opt, u);
    
    % Set up system matrix
    S = - fac_dt * T + C + A + B;
    S_II = S(grid.inner, grid.inner);
    S_IB = S(grid.inner, grid.DirichletBoundary);

    Tvold = -fac_dt * T * v_np1;
    rhs = rhs(grid.inner) - S_IB * v(grid.DirichletBoundary) + Tvold(grid.inner);

    % Solve system for inner vertices
    v_inner = S_II \ rhs;
    
    % Update v
    if strfind(prob.class, 'parabolic')
        v = sol.v;
    else
        % TODO: v.sol with correct boundary conditions is set for the parabolic case in solveParabolicProblem
        v = coef.g(grid.p');
    end
    v(grid.inner) = v_inner;

    % tmpsol.v = v;
    % tmpsol.u = u;
    % plotSolution(prob, grid, coef, opt, tmpsol);
    
    % Determine error in infinity norm
    error_v_inf = infError(v, v_prev)
    
    if error_v_inf < Howard_eps
        Howard_done = 1;
    else
        Howard_iter = Howard_iter + 1;
    end
    if Howard_iter > Howard_max_iter
        Howard_done = 1;
    end

end % while ~Howard_done

sol.v = v;
sol.u = u;

end % function [v] = solveNVPDE(prob, grid, coef, opt)
