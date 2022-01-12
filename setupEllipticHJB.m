function [p, c, g] = setupEllipticHJB(p, probLabel)

p.class = 'ellipticHJB';
p.DirichletBoundary = @(onBoundary,x) onBoundary;

if strcmpi(probLabel, 'MinimumArrivalTimex')
    p.dim = 2;
    p.mode = 'MinimumArrivalTimex';
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    p.control_dim = 1;

    % Setup grid
    g = setupGrid(p);

    % alpha = 1.5;
    alpha = 0;
    beta = 0.2;
    sigmax = 1;
    sigmay = 1;
    corrxy = 0.0;
    mux = 0.0;
    muy = 0.0;
    umax = 1;
    umin = -1;
    interest = 0.0;


	c.DirichletBC = 1;
    
	valx = 0.5*sigmax^2;
	valy = 0.5*sigmay^2;
	valxy = 0.5 * corrxy * sigmax * sigmay;
	
    % Determine length of input
    n = @(x) size(x,1);

	c.diffusion    = {@(x,u) valx*ones(n(x),1);
						@(x,u) valxy*ones(n(x),1);
						@(x,u) valxy*ones(n(x),1);
						@(x,u) valy*ones(n(x),1)};
	% c.diffusion    = {@(x) valx * x(:,1).^2 .* ones(n(x),1);
	% 									@(x) valxy*ones(n(x),1);
	% 									@(x) valxy*ones(n(x),1);
	% 									@(x) valy* x(:,2).^2 .* ones(n(x),1)};

	% c.diffusion    = {@(x) (0.5-x(:,1)).^2.*ones(n(x),1);
	%                   @(x) zeros(n(x),1);
	%                   @(x) zeros(n(x),1);
	%                   @(x) (0.5-x(:,2)).^2.*ones(n(x),1)};
	% c.diffusion    = {@(x) x(:,1).*ones(n(x),1);
	% 									@(x) zeros(n(x),1);
	% 									@(x) zeros(n(x),1);
	% 									@(x) (0.5-x(:,2)).^2.*ones(n(x),1)};
	
	c.convection   = {@(x,u) u(:,1) + mux *ones(n(x),1);
	                  @(x,u) muy *ones(n(x),1)};
	c.potential    = @(x,u) interest * ones(n(x),1);
    
	c.g  = @(x) zeros(n(x),1);
	c.f            = @(x,u) (1 + 0.5 * alpha * sqrt(u(:,1).^2) + beta * (abs(u(:,1))));
    if alpha > 0
        c.control_law = @(x, v, gradv, hessv) -1/alpha * ...
        [(gradv{1}<=-beta) .* (gradv{1}+beta) + (gradv{1}>beta) .* (gradv{1}-beta)];
    else
        c.control_law = @(x, v, gradv, hessv) ...
        [(gradv{1}<=-beta) * umax + (gradv{1}>beta) * umin];
    end
end % if strcmpi(probLabel, 'MinimumArrivalTimex')

if strcmpi(probLabel, 'MinimumArrivalTime')
    p.dim = 2;
    p.mode = 'MinimumArrivalTime';
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    p.control_dim = 2;

    % Setup grid
    g = setupGrid(p);

    alpha = 1.0;
    % alpha = 0;
    % beta = 0.2;
    beta = 0.0;
    sigmax = 1;
    sigmay = 1;
    corrxy = 0.0;
    mux = 0.0;
    muy = 0.0;
    umax = 1;
    umin = -1;
    interest = 0.0;


	c.DirichletBC = 1;
    
	valx = 0.5*sigmax^2;
	valy = 0.5*sigmay^2;
	valxy = 0.5 * corrxy * sigmax * sigmay;
	
    % Determine length of input
    n = @(x) size(x,1);

	c.diffusion    = {@(x,u) valx*ones(n(x),1);
						@(x,u) valxy*ones(n(x),1);
						@(x,u) valxy*ones(n(x),1);
						@(x,u) valy*ones(n(x),1)};
	% c.diffusion    = {@(x) valx * x(:,1).^2 .* ones(n(x),1);
	% 									@(x) valxy*ones(n(x),1);
	% 									@(x) valxy*ones(n(x),1);
	% 									@(x) valy* x(:,2).^2 .* ones(n(x),1)};

	% c.diffusion    = {@(x) (0.5-x(:,1)).^2.*ones(n(x),1);
	%                   @(x) zeros(n(x),1);
	%                   @(x) zeros(n(x),1);
	%                   @(x) (0.5-x(:,2)).^2.*ones(n(x),1)};
	% c.diffusion    = {@(x) x(:,1).*ones(n(x),1);
	% 									@(x) zeros(n(x),1);
	% 									@(x) zeros(n(x),1);
	% 									@(x) (0.5-x(:,2)).^2.*ones(n(x),1)};
	
	c.convection   = {@(x,u) u(:,1) + mux *ones(n(x),1);
	                  @(x,u) u(:,2) + muy *ones(n(x),1)};
	c.potential    = @(x,u) interest * ones(n(x),1);
    
	c.g  = @(x) zeros(n(x),1);
	c.f  = @(x,u) (1 + 0.5 * alpha * sqrt(u(:,1).^2 + u(:,2).^2) + beta * (abs(u(:,1)) + abs(u(:,2))));
    if alpha > 0
        c.control_law = @(x, v, gradv, hessv) -1/alpha * ...
        [(gradv{1}<=-beta) .* (gradv{1}+beta) + (gradv{1}>beta) .* (gradv{1}-beta), ...
         (gradv{2}<=-beta) .* (gradv{2}+beta) + (gradv{2}>beta) .* (gradv{2}-beta)];
    else
        c.control_law = @(x, v, gradv, hessv) ...
        [(gradv{1}<=-beta) * umax + (gradv{1}>beta) * umin, ...
         (gradv{2}<=-beta) * umax + (gradv{2}>beta) * umin];
    end
end % if strcmpi(probLabel, 'MinimumArrivalTime')

if strcmpi(probLabel, 'CinftyHJB')
    p.dim = 2;
    p.mode = 'CinftyHJB';
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    p.control_dim = 1;

    % Setup grid
    g = setupGrid(p);

	c.DirichletBC = 1;
    
    % Determine length of input
    n = @(x) size(x,1);

	c.diffusion    = {@(x,u) ones(n(x),1);
						@(x,u) zeros(n(x),1);
						@(x,u) zeros(n(x),1);
						@(x,u) ones(n(x),1)};
	
	c.convection   = {@(x,u) u(:,1);
	                  @(x,u) zeros(n(x),1)};

	c.potential    = @(x,u) zeros(n(x),1);

    c.solution = @(x) sin(x(:,1)*pi) .* sin(x(:,2)*pi);
    % c.optimal_control = @(x) x(:,1)-0.5;
    c.optimal_control = @(x) -3*exp(-x(:,1));
	
    c.g  = @(x) zeros(n(x),1);

    kappa = 0.4;
	c.f  = @(x,u) -1*(-2*pi^2*c.solution(x) ...
    - u .* pi .* (cos(pi*x(:,1)) .* sin(pi * x(:,2))) ...
    + kappa * (u - c.optimal_control(x)).^2);
    c.control_law = @(x,v,gradv,hessv) 1/(2*kappa) * (pi .* (cos(pi*x(:,1)) .* sin(pi * x(:,2))) - gradv{1}) + c.optimal_control(x);
end % if strcmpi(probLabel, 'CinftyHJB')

if strcmpi(probLabel, 'CinftyHJB2')
    p.dim = 2;
    p.mode = 'CinftyHJB2';
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    p.control_dim = 2;

    % Setup grid
    g = setupGrid(p);

	c.DirichletBC = 1;
    
    % Determine length of input
    n = @(x) size(x,1);

	c.diffusion    = {@(x,u) ones(n(x),1);
						@(x,u) zeros(n(x),1);
						@(x,u) zeros(n(x),1);
						@(x,u) ones(n(x),1)};
	

	
	c.convection   = {@(x,u) u(:,1);
	                  @(x,u) u(:,2)};

	c.potential    = @(x,u) zeros(n(x),1);

    c.sol = @(x) sin(x(:,1)*pi) .* sin(x(:,2)*pi);
    % c.optimal_control = @(x) x(:,1)-0.5;
    c.optimal_control = @(x) [-3*exp(-x(:,1)),sin(2*pi*x(:,2))];
    % c.optimal_control = @(x) zeros(n(x),2);
	
    c.g  = @(x) zeros(n(x),1);

    kappa = 4;
	c.f  = @(x,u) -1*(-2*pi^2*c.sol(x) ...
    - u(:,1) .* pi .* (cos(pi*x(:,1)) .* sin(pi * x(:,2))) ...
    - u(:,2) .* pi .* (sin(pi*x(:,1)) .* cos(pi * x(:,2))) ...
    + kappa * sum((u -c.optimal_control(x)).^2,2));
    % + kappa * (u(:,1) - (c.optimal_control(x))(:,1)).^2 ...
    % + kappa * (u(:,2) - (c.optimal_control(x))(:,2)).^2);
    c.control_law = @(x,v,gradv,hessv) ...
    [1/(2*kappa) * (pi .* (cos(pi*x(:,1)) .* sin(pi * x(:,2))) - gradv{1}), ...
     1/(2*kappa) * (pi .* (sin(pi*x(:,1)) .* cos(pi * x(:,2))) - gradv{2})] + c.optimal_control(x);

end % if strcmpi(probLabel, 'CinftyHJB')

% if strcmpi(probLabel, 'MertonInfiniteHorizon')
%     p.dim = 2;
%     p.mode = 'MertonInfiniteHorizon';
%     p.domain = 'UnitSquare';
%     [p.geom, p.gd] = getDomain(p.domain);
% 
%     % Set bounds of domain
%     p.xmin = 0.0;
%     p.xmax = 1.0;
%     p.ymin = 0.0;
%     p.ymax = 1.0;
% 
%     p.control_dim = 1;
% 
%     % Setup grid
%     g = setupGrid(p);
% 
%     % alpha = 1.5;
%     alpha = 0;
%     beta = 0.2;
%     
%     sigmax = 0.1;
%     sigmay = 0.2;
%     corrxy = 0.0;
% 
% 	valx = 0.5*sigmax^2;
% 	valy = 0.5*sigmay^2;
% 	valxy = 0.5 * corrxy * sigmax * sigmay;
% 
%     mux = 0.1;
%     muy = 0.13;
% 
%     interest = 0.02;
% 
% 	c.DirichletBC = 1;
%     
% 	
%     % Determine length of input
%     n = @(x) size(x,1);
% 
% 	c.diffusion    = {@(x,u) valx*(x(:,1).*u(:,1)).^2);
% 						@(x,u) valxy*(x(:,1).*u(:,1)).*(x(:,2).*u(:,2));
% 						@(x,u) valxy*(x(:,1).*u(:,1)).*(x(:,2).*u(:,2));
% 						@(x,u) valy*(x(:,2).*u(:,2)).^2};
% 	
% 	c.convection   = {@(x,u) (u(:,1)*(mux-interest) + interest)' *ones(n(x),1);
% 	                  @(x,u) muy *ones(n(x),1)};
% 	c.potential    = @(x,u) interest * ones(n(x),1);
%     
% 	c.g  = @(x) zeros(n(x),1);
% 	c.f            = @(x,u) (1 + 0.5 * alpha * sqrt(u(:,1).^2) + beta * (abs(u(:,1))));
%     if alpha > 0
%         c.control_law = @(x, v, gradv, hessv) -1/alpha * ...
%         [(gradv{1}<=-beta) .* (gradv{1}+beta) + (gradv{1}>beta) .* (gradv{1}-beta)];
%     else
%         c.control_law = @(x, v, gradv, hessv) ...
%         [(gradv{1}<=-beta) * umax + (gradv{1}>beta) * umin];
%     end
% end % if strcmpi(probLabel, 'MinimumArrivalTimex')

end % function [p, c, g] = setupEllipticHJB(p, probLabel)
