function [sol] = solveEllipticVI(prob, grid, coef, opt, sol)

if nargin < 5
    sol=[];
end
obstacle_problem = isa(coef.obstacle, 'function_handle');

% Check whether problem is parabolic or elliptic, this changes
% the behaviour of the code slightly.
if strfind(prob.class, 'parabolic')
    v_np1 = sol.v;
    fac_dt = 1/grid.dt;
    G = coef.g(grid.time, grid.p(:,grid.bnd)');

    if obstacle_problem
        Delta = coef.obstacle(grid.time, grid.p');
    end

    % if ~iscell(coef.obstacle)
    %     obs = coef.obstacle(grid.time, grid.p');
    % else
    %     m = length(coef.obstacle);
    % end
else
    v_np1 = zeros(grid.N,1);
    fac_dt = 0;
    G = coef.g(grid.p(:,grid.bnd)');
    
    if obstacle_problem
        Delta = coef.obstacle(grid.p');
    end


    % if ~iscell(coef.obstacle)
    %     obs = coef.obstacle(grid.p');
    % else
    %     m = length(coef.obstacle);
    % end
end


%% Set up parts of obstacle
if obstacle_problem
    Alpha = sparse(grid.N, grid.N);
    Beta = sparse(grid.N, grid.N);
    Gamma = speye(grid.N);

    % Initialize indicator of active set
    if isfield(sol,'pi')
        pi = sol.pi;
    else
        pi = zeros(grid.N,1);
    end
end

%% Set up parts of differential equation
% Get matrices for differential operators
A = getSecondOrderMatrix(prob, grid, coef);
B = getFirstOrderMatrix(prob, grid, coef, opt.upwinding);
C = getZerothOrderMatrix(prob, grid, coef);
T = grid.D0;
I = speye(grid.N);

% Differential operator
L = -fac_dt * T + A + B + C;
L(grid.bnd, :) = I(grid.bnd, :);

% Right-hand side
F = -fac_dt * T * v_np1 - setupLoadVector(grid, coef, opt);
F(grid.bnd) = G;

% Initialize iterator
converged = 0;
VI_EPS = 1e-13;
VI_iter = 1;
VI_iter_max = 1000;

% Start with zero value function
if isfield(sol,'v')
    v = sol.v;
else
    v = zeros(grid.N,1);
    v(grid.bnd) = G;
end
Dia = @(x) spdiags(x,[0],grid.N,grid.N);

% Compute residual for differential operator in non-divergence form

%% Non-trivial case
if ~obstacle_problem

    m = size(coef.obstacle.Alpha,1);
    t = grid.time;
    x = grid.p';
    for i=1:m
        Q{i} = getSecondOrderMatrix(prob, grid, coef.obstacle.Alpha{i}(t,x)) + ...
        getFirstOrderMatrix(prob, grid, coef.obstacle.Beta{i}(t,x), opt.upwinding) + ...
        getZerothOrderMatrix(prob, grid, coef.obstacle.Gamma{i}(t,x));

        Delta{i} = setupLoadVector(grid, coef.obstacle.Delta{i}(t,x), opt);
    end

end

fprintf(' Iter |   |res|_inf  |   |v-vold|_inf  \n')
fprintf('---------------------------------------\n')
if strcmpi(prob.obstacle_type, 'lower')
    Delta = -Delta;
    Gamma = -Gamma;
end

% Determine initial residuals
% R contains the residuals for the various parts in the max/min
% 1st component - 2nd order differential operator
R = L*v - F;
if obstacle_problem
    R = [R, (Alpha + Beta + Gamma)*v - Delta];
else
    
    % EPS = 1e-10;
    EPS = 0;
    for i=1:m
        R = [R, Q{i}*v - Delta{i}+i*EPS];
    end

end

% Determine maximum residual and maximizer
[res, pi] = prob.optimizer(R,[],2);


while ~converged

    vold = v;

    if obstacle_problem
        S = Dia(pi==1)*L + Dia(pi==2)*(Alpha + Beta + Gamma);
        rhs = Dia(pi==1)*F + Dia(pi==2)*Delta;
    else
        S = Dia(pi==1)*L;
        rhs = Dia(pi==1)*F;
        for i=1:m
            S = S + Dia(pi==i+1)*Q{i};
            rhs = rhs + Dia(pi==i+1)*Delta{i};
        end
    end

    % Enforce boundary conditions
    S(grid.bnd,:) = I(grid.bnd,:);
    rhs(grid.bnd) = G;

    v = S \ rhs;

    % Update residuals
    % R contains the residuals for the various parts in the max/min
    % 1st component - 2nd order differential operator
    R = L*v - F;
    if obstacle_problem
        R = [R, (Alpha + Beta + Gamma)*v - Delta];
    else
        
        % EPS = 1e-10;
        EPS = 0;
        for i=1:m
            R = [R, Q{i}*v - Delta{i}+i*EPS];
        end
    
    end

    % Determine maximum residual and maximizer
    [res, pi] = prob.optimizer(R,[],2);


    %% Regularization method %%
    % C =1e6;
    % 
    % Sreg = L + C * ( Dia(B1up*v>0)*B1up +  Dia(B2up*v>0)*B2up + ...
    %          Dia(S1up*v>0)*S1up +  Dia(S2up*v>0)*S2up);
    % rhsreg = F;
    % 
    % % Enforce boundary conditions
    % Sreg(grid.bnd,:) = L(grid.bnd,:);
    % rhsreg(grid.bnd) = F(grid.bnd);
    % 
    % v = Sreg\rhsreg;
    %% End regularization method %%
    

    fprintf(' %4d |  %6.4e  |  %6.4e \n', VI_iter, abs(prob.optimizer(res)), max(abs(v-vold)));
    if (abs(prob.optimizer(res)) < VI_EPS & VI_iter > 1) | VI_iter >= VI_iter_max
        converged = 1;
    end

    if max(abs(v-vold)) < VI_EPS & VI_iter > 1
        converged = 1;
    end

    if converged
        if strfind(prob.class, 'parabolic')
            if mod(100*grid.time,10) == 0 
                sfigure(2); clf; hold on
                plotOnTriangulation(grid.p,grid.t',v)
                sfigure(3); clf; hold on
                plotOnVoronoiGrid(grid, pi)
                keyboard
            end
        end
    end

    if ~converged
        if VI_iter > VI_iter_max
            converged = 1;
        else
            VI_iter = VI_iter + 1;
        end
    end
end

sol.v = v;
sol.pi = pi;

end % function [v] = solveNVPDE(prob, grid, coef, opt)
