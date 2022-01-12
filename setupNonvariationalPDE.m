function [p, c, g] = setupNonvariationalPDE(p, probLabel)

% Setup problem class specific data
p.class = 'ellipticPDE';
p.DirichletBoundary = @(onBoundary,x) onBoundary;

if strcmpi(probLabel,'ConvectionDominated')
    p.dim = 2;
    p.mode = 'ConvectionDominated';
    
    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    p.domain = 'UnitSquare';
    % p.domain = 'Diamond';
    [p.geom,p.gd] = getDomain(p.domain, p.xmin, p.xmax, p.ymin, p.ymax);


    % Setup grid
    g = setupGrid(p);

    mux = 0;
    muy = 0;

    % Determine length of input
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(x) 0.0*ones(n(x),1);

    % First order term
    c.convection   = {@(x) mux*ones(n(x),1);
    @(x) muy *ones(n(x),1)};
    
    corrxy = 0.0;
    % Second order term
    c.diffusion    = {@(x) ones(n(x),1);
    @(x) corrxy*ones(n(x),1);
    @(x) corrxy*ones(n(x),1);
    @(x) ones(n(x),1)};

    kappa = 1;
    c.solution = @(x) sin(kappa*pi*x(:,1)) .* sin(kappa*pi*x(:,2));
    c.ux = @(x) kappa*pi*cos(kappa*pi*x(:,1)) .* sin(kappa*pi*x(:,2));
    c.uy = @(x) kappa*pi*sin(kappa*pi*x(:,1)) .* cos(kappa*pi*x(:,2));
    c.uxx = @(x) -kappa^2*pi^2*sin(kappa*pi*x(:,1)) .* sin(kappa*pi*x(:,2));
    c.uxy = @(x) +kappa^2*pi^2*cos(kappa*pi*x(:,1)) .* cos(kappa*pi*x(:,2));
    % c.uyx = @(x) +kappa^2*pi^2*cos(kappa*pi*x(:,1)) .* cos(kappa*pi*x(:,2));
    c.uyy = @(x) -kappa^2*pi^2*sin(kappa*pi*x(:,1)) .* sin(kappa*pi*x(:,2));

    c.g = @(x) c.solution(x);
    c.f = @(x) -c.uxx(x) - c.uyy(x) - corrxy*(c.uxy(x)*2) - mux* c.ux(x) - muy*c.uy(x);
    % c.f = @(x) ones(n(x),1);
end

if strcmpi(probLabel,'Poisson2d')
    p.dim = 2;
    p.mode = 'Poisson';

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;
    p.domain = 'UnitSquare';
    [p.geom,p.gd] = getDomain(p.domain);
    % p.xmin = -1.0;
    % p.xmax = 1.0;
    % p.ymin = -1.0;
    % p.ymax = 1.0;
    % p.domain = 'Diamond';
    % [p.geom,p.gd] = getDomain(p.domain,p.xmin,p.xmax,p.ymin,p.ymax);

    % Setup grid
    g = setupGrid(p);


    % Formerly in setup coefficients
    n = @(x) size(x,1);
    R = 0.1;
    % Zeroth order term
    c.potential    = @(x) R*ones(n(x),1);
    mux = 1.;
    muy = -1.;
    % First order term
    c.convection   = {@(x) mux*ones(n(x),1);
    @(x) muy*ones(n(x),1)};

    % Second order term
    corrxy = +0.9;
    % Second order term
    c.diffusion    = {@(x) ones(n(x),1);
    @(x) corrxy*ones(n(x),1);
    @(x) corrxy*ones(n(x),1);
    @(x) ones(n(x),1)};

    c.g = @(x) zeros(n(x),1);
	% c.f = @(x) ones(n(x),1);
    % c.solution = @(x) 1/8*(x(:,1) .* (1-x(:,1)) +  x(:,2) .* (1-x(:,2)));
    c.solution = @(x) sin(2*pi*x(:,1)) .*sin(2*pi*x(:,2));
	c.f = @(x) (8*pi^2-R)*c.solution(x) - mux*2*pi*cos(2*pi*x(:,1)).*sin(2*pi*x(:,2)) - muy*2*pi*sin(2*pi*x(:,1)).*cos(2*pi*x(:,2)) - 8*pi^2*corrxy*cos(2*pi*x(:,1)).*cos(2*pi*x(:,2));
end

if strcmpi(probLabel,'Poisson1d')
    p.dim = 1;
    p.mode = 'Poisson';
    p.domain = 'UnitInterval';
    p.geom = [];

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;

    % Setup grid
    g = setupGrid(p);

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(x) zeros(n(x),1);

    % First order term
    c.convection   = {@(x) zeros(n(x),1)};

    % Second order term
    c.diffusion    = {@(x) ones(n(x),1)};

    c.g = @(x) zeros(n(x),1);
	c.f = @(x) +ones(n(x),1);
    c.solution = @(x) 0.5*x.*(1-x)
end

if strcmpi(probLabel,'SolInHalpha')
    p.alpha = 1.5;
    p.dim = 2;
    p.mode = 'SolInHalpha';
    p.domain = 'UnitSquare';
    [p.geom,p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    % Setup grid
    g = setupGrid(p);

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(x) zeros(n(x),1);

    % First order term
    c.convection   = {@(x) zeros(n(x),1); @(x) zeros(n(x),1)};

    % Second order term
    c.diffusion    = {@(x) ones(n(x),1); @(x) zeros(n(x),1); @(x) zeros(n(x),1); @(x) ones(n(x),1)};

    phi = @(x) atan2(x(:,1),x(:,2));
    c.solution = @(x) sqrt(x(:,1).^2 + x(:,2).^2).^(p.alpha) .* sin(2. * phi(x)) .* (1-x(:,1)) .* (1-x(:,2));

    vxx = @(x,y) (-1 + y) .* ((x.^2 + y.^2).^(p.alpha / 0.2e1 - 0.2e1)) ...
    .* (0.4e1 * ((x.^2 - x) * p.alpha + y.^2 + x) .* y .* cos(0.2e1 * atan2(x,y)) ...
    + sin(0.2e1 * atan2(x,y)) .* ((x.^ 2) .* (-1 + x) * p.alpha^2 + (x.^3 + 3 * x .* y.^2 + x.^2 - y.^2) ...
    .* p.alpha - (4 * y.^2 .* (-1 + x))));

    vyy = @(x,y) -0.4e1 * ((((-0.3e1 / 0.4e1 * p.alpha + 0.1e1) * y + p.alpha / 0.4e1 - 0.1e1) ...
    .* x.^2 - p.alpha * ((p.alpha + 0.1e1) .* y - p.alpha + 0.1e1) .* y.^2 / 0.4e1) ...
    .* sin(0.2e1 * atan2(x,y)) + (x.^2 + y .* (p.alpha * y - p.alpha + 0.1e1)) .* x ...
    .* cos(0.2e1 * atan2(x,y))) ...
    .* (-0.1e1 + x) .* (x.^2 + y.^2) .^ (p.alpha / 0.2e1 - 0.2e1);


    c.g = @(x) zeros(n(x),1);
	c.f = @(x) -vxx(x(:,1),x(:,2)) - vyy(x(:,1),x(:,2));

end

if strcmpi(probLabel,'CinftySol')
    p.dim = 2;
    p.mode = 'CinftySol';
    
    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    % p.domain = 'UnitSquare';
    p.domain = 'Diamond';
    
    [p.geom,p.gd] = getDomain(p.domain, p.xmin, p.xmax, p.ymin, p.ymax);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    % Setup grid
    g = setupGrid(p);

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(x) zeros(n(x),1);

    % First order term
    c.convection   = {@(x) zeros(n(x),1); @(x) zeros(n(x),1)};

    kappa = 0.3;
    % Second order term
    c.diffusion    = {@(x) ones(n(x),1); @(x) kappa*ones(n(x),1); @(x) kappa*ones(n(x),1); @(x) ones(n(x),1)};

    c.solution = @(x) sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2));

    c.g = @(x) c.solution(x);
	c.f = @(x) 8*pi^2*c.solution(x) - 2*kappa* 4*pi^2*cos(2*pi*x(:,1)) .* cos(2*pi*x(:,2));

end

if strcmpi(probLabel,'NonzeroBoundary')
    p.dim = 2;
    p.mode = 'NonzeroBoundary';
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    % Setup grid
    g = setupGrid(p);

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(x) zeros(n(x),1);

    % First order term
    c.convection   = {@(x) zeros(n(x),1); @(x) zeros(n(x),1)};

    kappa = 0.2
    % Second order term
    c.diffusion    = {@(x) ones(n(x),1); @(x) kappa*ones(n(x),1); @(x) kappa*ones(n(x),1); @(x) ones(n(x),1)};
    gamma = 1;
    c.solution = @(x) sin(gamma * x(:,1)) .* sin( gamma *x(:,2));

    c.g = @(x) c.solution(x);
	c.f = @(x) 2 * gamma^2 *c.solution(x) - 2*kappa* gamma^2*cos(gamma*x(:,1)) .* cos(gamma*x(:,2));

end

if strcmpi(probLabel,'LognormalDist')
    p.dim = 2;
    p.mode = 'LognormalDist';
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    mux = 0.8;
    muy = 0.5;
    sigmax = 0.4;
    sigmay = 0.62;
    corrxy = 0.2;

    % Setup grid
    g = setupGrid(p);

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(x) 0.1*ones(n(x),1);

    % First order term
    c.convection   = {@(x) mux * x(:,1); @(x) muy * x(:,2)};
    % c.convection   = {@(x) mux * ones(n(x),1); @(x) muy * ones(n(x),1)};

	valx = 0.5*sigmax^2;
	valy = 0.5*sigmay^2;
	valxy = 0.5 * corrxy * sigmax * sigmay;
    % Second order term
    c.diffusion    = {@(x) valx *x(:,1).^2; ...
                      @(x) valxy*x(:,1).*x(:,2); ...
                      @(x) valxy*x(:,1).*x(:,2); ...
                      @(x) valy *x(:,2).^2};
    % c.diffusion    = {@(x) 0.01+0.5*(sqrt(x(:,1).^2 + x(:,2).^2)<0.6); ...
    %                   @(x) zeros(n(x),1); ...
    %                   @(x) zeros(n(x),1); ...
    %                   @(x) 0.01+0.5*(sqrt(x(:,1).^2 + x(:,2).^2)<0.6)};

    gamma = pi;
    
    c.solution = @(x) sin(gamma * x(:,1)) .* sin( gamma *x(:,2));
    
    c.vx = @(x) gamma * cos(gamma * x(:,1)) .* sin( gamma *x(:,2));
    c.vy = @(x) gamma * sin(gamma * x(:,1)) .* cos( gamma *x(:,2));

    c.vxx = @(x) -gamma^2 * sin(gamma * x(:,1)) .* sin( gamma *x(:,2));
    c.vxy = @(x) gamma^2 * cos(gamma * x(:,1)) .* cos( gamma *x(:,2));
    c.vyy = @(x) -gamma^2 * sin(gamma * x(:,1)) .* sin( gamma *x(:,2));

    c.g = @(x) c.solution(x);
	c.f = @(x) -1*(c.diffusion{1}(x) .* c.vxx(x) + ...
              (c.diffusion{2}(x) + c.diffusion{3}(x)) .* c.vxy(x) + ...
               c.diffusion{4}(x) .* c.vyy(x) + ...
               c.convection{1}(x) .* c.vx(x) + ...
               c.convection{2}(x) .* c.vy(x) + ...
               c.potential(x) .* c.solution(x));

end

end % function [p, c, g] = setupNonvariationalPDE(p, probLabel)
