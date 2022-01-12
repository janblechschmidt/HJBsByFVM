function [sol] = solveNVPDE(prob, grid, coef, opt, sol)

% Initialize solution vector with correct boundary conditions
if strfind(prob.class, 'parabolic')
    v = coef.g(grid.time, grid.p');
    vold = sol.v;
    fac_dt = 1/grid.dt;
else
    v = coef.g(grid.p');
    vold = zeros(size(v));
    fac_dt = 0;
end

% Assemble coeffcient matrices
C = getZerothOrderMatrix(prob, grid, coef, []);
T = grid.D0;
B = getFirstOrderMatrix(prob, grid, coef, opt.upwinding, []);
A = getSecondOrderMatrix(prob, grid, coef, []);

% Assemble right-hand side
rhs = - setupLoadVector(grid, coef, opt);
% Set up system matrix
S = - fac_dt * T + C + A + B;
% Sback = S;
% keyboard
S_II = S(grid.inner, grid.inner);
% if length(grid.OutflowBoundary) > 0
%     bcoef = evaluateCoefficient(prob, grid, coef.convection, []);
%     rhs_outflow = 0;
%     for i=1:prob.dim
%         rhs_outflow = rhs_outflow + grid.Dout{i}*bcoef{i};
%     end
% end
% Check M-matrix property of S_II

if isfield(prob,'Tmax') & grid.time == .999
    diagvals = abs(diag(S_II));
    S_IIdiagZero = S_II - diag(diag(S_II));
    offdiagvals = sum(abs(S_IIdiagZero),2);
    fprintf('Test diagonal dominance of S_II: %6.4e\n',full(min(diagvals-offdiagvals)));
    % pause
end

S_IB = S(grid.inner, grid.DirichletBoundary);
% S_IB = S(grid.InnerIndices, grid.bnd);
Tvold = -fac_dt * T * vold;
rhs = rhs(grid.inner) - S_IB * v(grid.DirichletBoundary) + Tvold(grid.inner);
% rhs = rhs(grid.InnerIndices) - S_IB * v(grid.bnd) + Tvold(grid.InnerIndices);

% Solve system for inner vertices
v_inner = S_II \ rhs;

% Incorporate values at inner nodes into v
v(grid.inner) = v_inner;

sol.v = v;

if opt.plotGradient

    gradv = cellfun(@(C) C*v,grid.D1,'UniformOutput',false');
    % Plot gradient
    sfigure(7); clf;
    for i = 1:prob.dim
        subplot(prob.dim, 1, i);
        plotOnTriangulation(grid.p, grid.t' ,gradv{i});
        title('gradv')
    end
end
if opt.plotHessian

    hessv = cellfun(@(C) C*v,grid.D2,'UniformOutput',false');
    % Plot Hessian matrix
    sfigure(8); clf;
    for i = 1:prob.dim
        for j = 1:prob.dim
            subplot(prob.dim, prob.dim, 2*(i-1)+j);
            plotOnTriangulation(grid.p, grid.t' ,gradv{i});
            title(sprintf('D^2_{x%ix%i} v', i, j));
        end
    end
end
if opt.plotGradient | opt.plotHessian
    fprintf('Hit a key to proceed...\n');
    pause
end
end % function [sol] = solveNVPDE(prob, grid, coef, opt, dt, vold)
