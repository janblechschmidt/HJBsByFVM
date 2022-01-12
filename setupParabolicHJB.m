function [p, c, g] = setupParabolicHJB(p, probLabel)

% Setup problem class specific data
p.class = 'parabolicHJB';
p.DirichletBoundary = @(onBoundary,x) onBoundary;

if strcmpi(probLabel, 'ParabolicMinimumArrivalTime')
    p.dim = 2;
    p.mode = 'ParabolicMinimumArrivalTime';
    % p.domain = 'UnitSquare';
    p.domain = 'Rectangle';

    % Set bounds of domain
    p.xmin = -1.0;
    p.xmax = 1.0;
    p.ymin = -1.0;
    p.ymax = 1.0;
    p.Tmin = 0.0;
    p.Tmax = 1.0;
    [p.geom, p.gd] = getDomain(p.domain, p.xmin, p.xmax, p.ymin, p.ymax);

    p.control_dim = 2;

    % Setup grid
    g = setupGrid(p);

    if ~isfield(p, 'alpha')
        p.alpha = 0.5;
    end

    if ~isfield(p, 'beta')
        p.beta = 0.5;
    end
    % p.alpha = 1.0;
    % p.beta = 0.0;
    % p.alpha = 0;
    % p.beta = 0.2;
    mux = 0.0;
    muy = 0.0;
    p.Umax = 1;
    p.Umin = -1;
    discount = 0.0;


	c.DirichletBC = 1;
    
    if ~isfield(p, 'sigmax')
        p.sigmax = 0.5;
    end
    if ~isfield(p, 'sigmay')
        p.sigmay = 0.5;
    end
    if ~isfield(p, 'corrxy')
        p.corrxy = 0.0;
    end

	valx = 0.5*p.sigmax^2;
	valy = 0.5*p.sigmay^2;
	valxy = 0.5 * p.corrxy * p.sigmax * p.sigmay;
	
    % Determine length of input
    n = @(x) size(x,1);

	% c.diffusion    = {@(t, x,u) valx*ones(n(x),1);
	%                     @(t, x,u) valxy*ones(n(x),1);
	%                     @(t, x,u) valxy*ones(n(x),1);
	%                     @(t, x,u) valy*ones(n(x),1)};
    c.diffusion    = {@(t, x,u) valx*ones(n(x),1);
                        @(t, x, u) (2*(sqrt(x(:,1).^2+x(:,2).^2)<.5)-1).*valxy.*ones(n(x),1);
                        @(t, x, u) (2*(sqrt(x(:,1).^2+x(:,2).^2)<.5)-1).*valxy.*ones(n(x),1);
                        @(t, x,u) valy*ones(n(x),1)};

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
	
	c.convection   = {@(t, x, u) u(:,1) + mux *ones(n(x),1);
	                  @(t, x, u) u(:,2) + muy *ones(n(x),1)};
	c.potential    = @(t, x, u) discount * ones(n(x),1);
    
	c.g  = @(t, x) zeros(n(x),1);
	c.f  = @(t, x,u) 1 + 0.5 * p.alpha * sum(u.^2, 2) + p.beta * sum(abs(u),2);

    c.control_law = @control_MinimumArrivalTime;
    % if p.alpha > 0
    %     c.control_law = @(t, x, v, gradv, hessv) -1/p.alpha * ...
    %     [(gradv{1}<=-p.beta) .* (gradv{1}+p.beta) + (gradv{1}>p.beta) .* (gradv{1}-p.beta), ...
    %      (gradv{2}<=-p.beta) .* (gradv{2}+p.beta) + (gradv{2}>p.beta) .* (gradv{2}-p.beta)];
    % else
    %     c.control_law = @(t, x, v, gradv, hessv) ...
    %     [(gradv{1}<=-p.beta) * umax + (gradv{1}>p.beta) * p.Umin, ...
    %      (gradv{2}<=-p.beta) * umax + (gradv{2}>p.beta) * p.Umin];
    % end
end

if strcmpi(probLabel,'OptimalPortfolio2d')

    % Transaction cost
    gamma = 0.5;
    % Utility exponent
    alpha = 0.3;
    % Interest rate
    r = 0.07;
    % Drift
    mu = [0.15;0.02];
    % Diffusion
    sigma = [0.42, 0.1; 0.1, 0.38];
    A = sigma * sigma';
    
    p.dim = 2;
    % p.mode = 'Poisson';
    p.control_dim = 1;

    % Spatial dimensions
    p.xmin = -1/gamma;
    p.xmax = 1/gamma;
    p.ymin = -1/gamma;
    p.ymax = 1/gamma;

    % Time horizon
    p.Tmax = 1.0;
    p.Tmin = 0.0;

    p.domain = 'Diamond';
    prob.mesh_type = 'irregular';
    [p.geom,p.gd] = getDomain(p.domain,p.xmin,p.xmax,p.ymin,p.ymax);
    if ~isfield(p, 'hmax')
        p.hmax = 0.1;
    end
    if ~isfield(p, 'dtmax')
        p.dtmax = 0.1;
    end

    % Setup grid
    g = setupGrid(p);

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(t,x,u) alpha * (r + x*(mu-r) + ...
        0.5*(1-alpha)*sum(sum((x*A).*x,2)));
    
    % First order term
    c.convection   = {@(t,x,u) ( mu(1)-r - ( (mu(1)-r) * x(:,1) + (mu(2) - r) * x(:,2)  )  ) .*x(:,1) ...
        -  0.5 * (1-alpha) * ( A(1,1)*x(:,1).^2      .* (2-2*x(:,1)) ...
                           + 2*A(1,2)*x(:,1).*x(:,2) .* (1-2*x(:,1)) ...
                           +   A(2,2)*x(:,2).^2      .* (0-2*x(:,1)));
    @(t,x,u) ( mu(2)-r - ( (mu(1)-r) * x(:,1) + (mu(2) - r) * x(:,2)  )  ) .*x(:,2) ...
        -  0.5 * (1-alpha) * ( A(1,1)*x(:,1).^2      .* (0-2*x(:,2)) ...
                           + 2*A(1,2)*x(:,1).*x(:,2) .* (1-2*x(:,2)) ...
                           +   A(2,2)*x(:,2).^2      .* (2-2*x(:,2))) };

    % Second order term
    c.diffusion    = {@(t,x,u) -0.5*( A(1,1) * x(:,1).^2      .* (1-x(:,1)).^2 ...
                                  + A(1,2) * x(:,1).*x(:,2) .* (-x(:,1).*(1-x(:,1)) ) ...
                                  + A(2,1) * x(:,2).*x(:,2) .* (-x(:,1).*(1-x(:,1)) ) ...
                                  + A(2,2) * x(:,2).^2      .* (-x(:,1)).^2);
                      @(t,x,u) -0.5*( A(1,1) * x(:,1).^2      .* (1-x(:,1)).*(-x(:,2)) ...
                                  + A(1,2) * x(:,1).*x(:,2) .* (1-x(:,1)).*(1-x(:,2)) ...
                                  + A(2,1) * x(:,2).*x(:,2) .* (-x(:,1)) .*(-x(:,2)) ...
                                  + A(2,2) * x(:,2).^2      .* (-x(:,1)) .*(1-x(:,2)) );
                      @(t,x,u) -0.5*( A(1,1) * x(:,1).^2      .* (1-x(:,1)).*(-x(:,2)) ...
                                  + A(1,2) * x(:,1).*x(:,2) .* (-x(:,1)).*(-x(:,2)) ...
                                  + A(2,1) * x(:,2).*x(:,2) .* (1-x(:,1)) .*(1-x(:,2)) ...
                                  + A(2,2) * x(:,2).^2      .* (-x(:,1)) .*(1-x(:,2)) );
                      @(t,x,u) -0.5*( A(1,1) * x(:,1).^2      .* (-x(:,2)).^2 ...
                                  + A(1,2) * x(:,1).*x(:,2) .* (-x(:,2).*(1-x(:,2)) ) ...
                                  + A(2,1) * x(:,2).*x(:,2) .* (-x(:,2).*(1-x(:,2)) ) ...
                                  + A(2,2) * x(:,2).^2      .* (1-x(:,2)).^2)};

    c.control_law = @(t, x, v, gradv, hessv)  zeros(n(x),1);
    c.g = @(t,x) 1/alpha * ( 1 - gamma * (abs(x(:,1)) + abs(x(:,2))) ).^alpha;
	c.f = @(t,x,u) zeros(n(x),1);
end

if strcmpi(probLabel, 'MATCircle')
    p.dim = 2;
    p.mode = 'ParabolicMinimumArrivalTime';
    load('Circdom.mat')
    p.geom = g;
    p.gd = gd;
    % p.domain = 'UnitSquare';
    % p.domain = 'Rectangle';
    % Set bounds of domain
    % p.xmin = -1.0;
    % p.xmax = 1.0;
    % p.ymin = -1.0;
    % p.ymax = 1.0;
    p.Tmin = 0.0;
    p.Tmax = 1.0;
    % [p.geom, p.gd] = getDomain(p.domain, p.xmin, p.xmax, p.ymin, p.ymax);

    p.control_dim = 2;

    % Setup grid
    g = setupGrid(p,pnts,e,t);

    if ~isfield(p, 'alpha')
        p.alpha = 0.5;
    end

    if ~isfield(p, 'beta')
        p.beta = 0.5;
    end
    % p.alpha = 1.0;
    % p.beta = 0.0;
    % p.alpha = 0;
    % p.beta = 0.2;
    mux = 0.0;
    muy = 0.0;
    p.Umax = 1;
    p.Umin = -1;
    discount = 0.0;


	c.DirichletBC = 1;
    
    if ~isfield(p, 'sigmax')
        p.sigmax = 0.5;
    end
    if ~isfield(p, 'sigmay')
        p.sigmay = 0.5;
    end
    if ~isfield(p, 'corrxy')
        p.corrxy = 0.0;
    end
	valx = 0.5*p.sigmax^2;
	valy = 0.5*p.sigmay^2;
	valxy = 0.5 * p.corrxy * p.sigmax * p.sigmay;
	
    % Determine length of input
    n = @(x) size(x,1);

	c.diffusion    = {@(t, x,u) valx*ones(n(x),1);
						@(t, x,u) valxy*ones(n(x),1);
						@(t, x,u) valxy*ones(n(x),1);
						@(t, x,u) valy*ones(n(x),1)};
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
	
	c.convection   = {@(t, x, u) u(:,1) + mux *ones(n(x),1);
	                  @(t, x, u) u(:,2) + muy *ones(n(x),1)};
	c.potential    = @(t, x, u) discount * ones(n(x),1);
    
	c.g  = @(t, x) zeros(n(x),1);
	c.f  = @(t, x,u) 1 + 0.5 * p.alpha * sum(u.^2, 2) + p.beta * sum(abs(u),2);

    c.control_law = @control_MinimumArrivalTime;
    % if p.alpha > 0
    %     c.control_law = @(t, x, v, gradv, hessv) -1/p.alpha * ...
    %     [(gradv{1}<=-p.beta) .* (gradv{1}+p.beta) + (gradv{1}>p.beta) .* (gradv{1}-p.beta), ...
    %      (gradv{2}<=-p.beta) .* (gradv{2}+p.beta) + (gradv{2}>p.beta) .* (gradv{2}-p.beta)];
    % else
    %     c.control_law = @(t, x, v, gradv, hessv) ...
    %     [(gradv{1}<=-p.beta) * umax + (gradv{1}>p.beta) * p.Umin, ...
    %      (gradv{2}<=-p.beta) * umax + (gradv{2}>p.beta) * p.Umin];
    % end
end
