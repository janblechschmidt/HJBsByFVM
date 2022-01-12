function [p, c, g] = setupNonvariationalPDE(p, probLabel)

% Setup problem class specific data
n = @(x) size(x,1);
p.class = 'ellipticVI';
p.mode = probLabel;
p.DirichletBoundary = @(onBoundary,x) onBoundary;
p.optimizer = @max;
p.obstacle_type = 'upper';

if strcmpi(probLabel,'PoissonObstacle')
    p.dim = 2;
    p.domain = 'UnitSquare';
    [p.geom, p.gd] = getDomain(p.domain);

    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 1.0;
    p.ymin = 0.0;
    p.ymax = 1.0;

    % Setup grid
    g = setupGrid(p);

    % Zeroth order term
    c.potential    = @(x) zeros(n(x),1);

    % First order term
    c.convection   = {@(x) zeros(n(x),1);
    @(x) zeros(n(x),1)};

    % Second order term
    c.diffusion    = {@(x) -ones(n(x),1);
    @(x) zeros(n(x),1);
    @(x) zeros(n(x),1);
    @(x) -ones(n(x),1)};

    c.obstacle = @(x) 0.04*ones(n(x),1);
    c.g = @(x) zeros(n(x),1);
	c.f = @(x) -ones(n(x),1);
    % c.solution = @(x) 1/8*(x(:,1) .* (1-x(:,1)) +  x(:,2) .* (1-x(:,2)));
end


end % function [p, c, g] = setupNonvariationalPDE(p, probLabel)
