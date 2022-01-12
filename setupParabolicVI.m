function [p, c, g] = setupParabolicVI(p, probLabel)

% Short-hand notation for size of x
n = @(x) size(x,1);

p.class = 'parabolicVI'
p.DirichletBoundary = @(onBoundary,x) onBoundary;

p.mode = probLabel;
p.optimizer = @max;
p.obstacle_type = 'upper';

if strcmpi(probLabel,'AmericanWorstOfTwoAssetsPut')
    p.dim = 2;
    p.domain = 'Rectangle';
    % Set bounds of domain
    p.xmin = 0.0;
    p.xmax = 80.0;
    p.ymin = 0.0;
    p.ymax = 80.0;

    % The geometry is not necessary for the regular case. I should try to avoid this anyways.
    [p.geom, p.gd] = getDomain(p.domain, p.xmin, p.xmax, p.ymin, p.ymax);

    % Dirichlet boundary
    p.DirichletBoundary = @(onBoundary,x) onBoundary & x(:,1) >= p.xmax-1e-5 & x(:,2) >= p.ymax - 1e-5;

    p.obstacle_type = 'lower';

    % Overwrite hmax
    % p.hmax = 0.2;
    % p.dtmax = 0.025;
    % p.hmax = 0.1;
    % p.dtmax = 0.01;

    p.hmax = 0.5;
    p.dtmax = 0.01;
    % p.hmax = 10;
    % p.dtmax = 0.1;

    % Time horizon
    p.Tmax = 1.0;
    p.Tmin = 0.5;

    r = 0.05;
    sigma1 = 0.3;
    sigma2 = 0.3;
    rho = 0.5;
    K = 40;

    % Zeroth order term
    mode = 'zvan'
    switch mode
    case 'simple'

        c.potential    = @(t,x) -r * ones(n(x),1);
        c.convection   = {@(t,x) +r * x(:,1);
                          @(t,x) +r * x(:,2)};
        c.diffusion    = {@(t,x) +0.5*sigma1^2 * x(:,1).^2;
                          @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2);
                          @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2);
                          @(t,x) +0.5*sigma2^2 * x(:,2).^2};
    case 'own'
        c.potential    = @(t,x) -r * ones(n(x),1);
        c.convection   = {@(t,x) +r * x(:,1);
                          @(t,x) +r * x(:,2)};
        c.diffusion    = {@(t,x) +0.5*sigma1^2 * x(:,1).^2;
                          @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2);
                          @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2);
                          @(t,x) +0.5*sigma2^2 * x(:,2).^2};
        % Dirichlet boundary
        p.DirichletBoundary = @(onBoundary,x) onBoundary & ...
            ( (x(:,1) >= p.xmax-1e-5 & x(:,2) >= p.ymax - 1e-5) | ...
            (x(:,1) < p.xmin+1e-5) | (x(:,2) < p.ymin+1e-5) );

    case 'zvan'
        c.potential    = @(t,x) -r * ones(n(x),1);
        % First order term
        c.convection   = {@(t,x) +r * x(:,1) ...
                          .* (x(:,1) > p.xmin) .* (x(:,1)<p.xmax);
                          @(t,x) +r * x(:,2) ...
                          .* (x(:,2) > p.ymin) .* (x(:,2) < p.ymax)};

        % Second order term
        c.diffusion    = {@(t,x) +0.5*sigma1^2 * x(:,1).^2 ...
                       .* (x(:,1) > p.xmin) .* (x(:,1)<p.xmax);
                       @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2) ...
                       .* (x(:,1) > p.xmin) .* (x(:,1)<p.xmax) .* (x(:,2) > p.ymin) .* (x(:,2) < p.ymax);
                       @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2) ...
                       .* (x(:,1) > p.xmin) .* (x(:,1)<p.xmax) .* (x(:,2) > p.ymin) .* (x(:,2) < p.ymax);
                       @(t,x) +0.5*sigma2^2 * x(:,2).^2 ...
                       .* (x(:,2) > p.ymin) .* (x(:,2) < p.ymax)};
    end

    c.obstacle = @(t,x) max(K-min(x(:,1),x(:,2)),0);
    c.g = @(t,x) exp(-r*(p.Tmax-t))*max(K-min(x(:,1),x(:,2)),0);
	c.f = @(t,x) zeros(n(x),1);
    
end

if strcmpi(probLabel,'OptimalPortfolio2d')

    % Transaction cost
    gamma = 0.025;
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
    p.mode = 'Poisson';
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
    [p.geom,p.gd] = getDomain(p.domain,p.xmin,p.xmax,p.ymin,p.ymax);

    p.hmax = 2;
    p.dtmax = 0.01;

    % Formerly in setup coefficients
    n = @(x) size(x,1);

    % Zeroth order term
    c.potential    = @(t,x) alpha * (r + x*(mu-r) - ...
        0.5*(1-alpha)*sum(sum((x*A).*x,2)));
    
    % First order term
    c.convection   = {
    @(t,x) ( mu(1)-r - ( (mu(1)-r) * x(:,1) + (mu(2) - r) * x(:,2)  )  ) .*x(:,1) ...
        -  0.5 * (1-alpha) * ( A(1,1)*x(:,1).^2             .* (2-2*x(:,1)) ...
                           + (A(1,2)+A(2,1))*x(:,1).*x(:,2) .* (1-2*x(:,1)) ...
                           +   A(2,2)*x(:,2).^2             .* (0-2*x(:,1)));
    @(t,x) ( mu(2)-r - ( (mu(1)-r) * x(:,1) + (mu(2) - r) * x(:,2)  )  ) .*x(:,2) ...
        -  0.5 * (1-alpha) * ( A(1,1)*x(:,1).^2             .* (0-2*x(:,2)) ...
                           + (A(1,2)+A(2,1))*x(:,1).*x(:,2) .* (1-2*x(:,2)) ...
                           +   A(2,2)*x(:,2).^2             .* (2-2*x(:,2))) };

    % Second order term
    c.diffusion    = {@(t,x) +0.5*( A(1,1) * x(:,1).^2      .* (1-x(:,1)).^2 ...
                                  - (A(1,2) + A(2,1)) * x(:,1).^2 .* (1-x(:,1)) .* x(:,2)  ...
                                  + A(2,2) * x(:,2).^2 .* x(:,1).^2);
                      @(t,x) +0.5*( -A(1,1) * x(:,1).^2 .* (1-x(:,1)) .*(x(:,2)) ...
                                  + A(1,2) * x(:,1).*x(:,2) .* (1-x(:,1)).*(1-x(:,2)) ...
                                  + A(2,1) * x(:,2).^2.*x(:,1).^2 ...
                                  - A(2,2) * x(:,2).^2 .* x(:,1) .*(1-x(:,2)) );
                      @(t,x) +0.5*( -A(1,1) * x(:,1).^2 .* (1-x(:,1)) .*(x(:,2)) ...
                                  + A(1,2) * x(:,2).^2.*x(:,1).^2 ...
                                  + A(2,1) * x(:,1).*x(:,2) .* (1-x(:,1)).*(1-x(:,2)) ...
                                  - A(2,2) * x(:,2).^2 .* x(:,1) .*(1-x(:,2)) );
                      @(t,x) +0.5*( A(1,1) * x(:,1).^2      .* x(:,2).^2 ...
                                  - (A(1,2)+A(2,1)) * x(:,1).*x(:,2).^2 .*(1-x(:,2)) ...
                                  + A(2,2) * x(:,2).^2      .* (1-x(:,2)).^2)};

    c.obstacle.Alpha = {@(t,x) zeros(n(x),4); @(t,x) zeros(n(x),4); @(t,x) zeros(n(x),4); @(t,x) zeros(n(x),4) };
    c.obstacle.Beta = {@(t,x) [1 + gamma * x(:,1), gamma * x(:,2)];
                    @(t,x) [gamma * x(:,1), 1 + gamma * x(:,2)];
                    @(t,x) - [1 - gamma * x(:,1), -gamma * x(:,2)];
                    @(t,x) - [-gamma * x(:,1), 1 - gamma * x(:,2)]};
    c.obstacle.Gamma = {@(t,x) -alpha*gamma*ones(n(x),1);
                    @(t,x) -alpha*gamma*ones(n(x),1);
                    @(t,x) -alpha*gamma*ones(n(x),1);
                    @(t,x) -alpha*gamma*ones(n(x),1)};
    c.obstacle.Delta = {@(t,x) zeros(n(x),1); @(t,x) zeros(n(x),1); @(t,x) zeros(n(x),1); @(t,x) zeros(n(x),1) };

    c.g = @(t,x) 1/alpha * ( 1 - gamma * (abs(x(:,1)) + abs(x(:,2))) ).^alpha;
	c.f = @(t,x) zeros(n(x),1);
end

% Setup grid
g = setupGrid(p);

end % function [p, c, g] = setupNonvariationalPDE(p, probLabel)
