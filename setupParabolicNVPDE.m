function [p, c, g] = setupParabolicNVPDE(p, probLabel)

    % Short-hand notation for size of x
    n = @(x) size(x,1);

    p.DirichletBoundary = @(onBoundary,x) onBoundary;

    p.mode = probLabel;
    p.class = 'parabolicPDE'
    if strcmpi(probLabel,'CinftySolTime')
        p.dim = 2;
        p.domain = 'UnitSquare';
        [p.geom, p.gd] = getDomain(p.domain);

        % Set bounds of domain
        p.xmin = 0.0;
        p.xmax = 1.0;
        p.ymin = 0.0;
        p.ymax = 1.0;

        % Time horizon
        p.Tmax = 1.0;
        p.Tmin = 0.0;


        % Zeroth order term
        c.potential    = @(t,x) zeros(n(x),1);

        % First order term
        c.convection   = {@(t,x) zeros(n(x),1); @(t,x) zeros(n(x),1)};

        kappa = 0.3;
        % Second order term
        c.diffusion    = {@(t,x) ones(n(x),1); @(t,x) kappa*ones(n(x),1); @(t,x) kappa*ones(n(x),1); @(t,x) ones(n(x),1)};

        gamma = -0.5;
        c.solution = @(t,x) exp(gamma * t) * sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2));

        c.g = @(t,x) zeros(n(x),1);
        c.f = @(t,x) (8*pi^2 - gamma) *c.solution(t,x) - 2*kappa* 4 * pi^2 * exp(gamma * t) * cos(2*pi*x(:,1)) .* cos(2*pi*x(:,2));

    end

    if strcmpi(probLabel,'WorstOfTwoAssetsPut')
        p.dim = 2;

        % Set bounds of domain
        p.xmin = 0.0;
        p.xmax = 140.0;
        p.ymin = 0.0;
        p.ymax = 140.0;

        p.domain = 'Rectangle';
        % The geometry is not necessary for the regular case. I should try to avoid this anyways.
        [p.geom, p.gd] = getDomain(p.domain, p.xmin, p.xmax, p.ymin, p.ymax);

        % Dirichlet boundary
        p.DirichletBoundary = @(onBoundary,x) onBoundary & x(:,1) >= p.xmax-1e-9 & x(:,2) >= p.ymax - 1e-9;

        % Time horizon
        p.Tmax = 1.0;
        p.Tmin = 0.5;

        r = 0.05;
        sigma1 = 0.3;
        sigma2 = 0.3;
        rho = 0.5;
        K = 40;

        % Zeroth order term
        mode = 'own';
        % mode = 'zvan';
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
                % c.convection   = {@(t,x) +r * x(:,1) .* (x(:,1) < p.xmax);
                %     @(t,x) +r * x(:,2) .* (x(:,2) < p.ymax)};
                c.diffusion    = {@(t,x) +0.5*sigma1^2 * x(:,1).^2 .* (x(:,1)<p.xmax);
                    @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2) .* (x(:,1)<p.xmax).* (x(:,2)<p.ymax);
                    @(t,x) +0.5*rho*sigma1*sigma2*x(:,1).*x(:,2) .* (x(:,1)<p.xmax).* (x(:,2)<p.ymax);
                    @(t,x) +0.5*sigma2^2 * x(:,2).^2 .* (x(:,2) < p.ymax)};
                % Dirichlet boundary
                p.DirichletBoundary = @(onBoundary,x) onBoundary & ...
                    ( (x(:,1) >= p.xmax-1e-9 & x(:,2) >= p.ymax - 1e-9) | ...
                    (x(:,1) < p.xmin+1e-9) | (x(:,2) < p.ymin+1e-9) );

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

        c.g = @(t,x) exp(-r*(p.Tmax-t))*max(K-min(x(:,1),x(:,2)),0);
        c.f = @(t,x) zeros(n(x),1);

        c.solution = explicitSolutionWorstOfTwoAssetsPut(K,r,sigma1,sigma2,rho,p.Tmax);
    end

    if strcmpi(probLabel,'AsianCall')
        p.dim = 2;
        p.domain = 'UnitSquare';
        [p.geom, p.gd] = getDomain(p.domain);

        % Set bounds of domain
        p.xmin = 0.0;
        p.xmax = 200.0;
        p.ymin = 0.0;
        p.ymax = 200.0;

        % Overwrite initial setting of h and dt
        % p.hmax = 4;
        % p.dtmax = 0.02;
        p.hmax = 2;
        p.dtmax = 0.01;
        % p.hmax = 1;
        % p.dtmax = 0.005;

        % Time horizon
        %p.Tmax = 1.0;
        %p.Tmin = 0.75;
        p.Tmax = 0.25;
        p.Tmin = 0.0;

        % Dirichlet boundary
        p.DirichletBoundary = @(onBoundary,x) onBoundary & x(:,1) >= p.xmax-1e-9 & x(:,2) >= p.ymax - 1e-9;

        % Values from Zvan
        sigmaS = 0.1;
        r = 0.1;
        K = 100;

        % Zeroth order term
        c.potential    = @(t,x) -r * ones(n(x),1);
        % .* (x(:,1)<p.xmax);

        % First order term
        c.convection   = {@(t,x) r*x(:,1) .* (x(:,1) > p.xmin) .* (x(:,1)<p.xmax);
            @(t,x) posConvection(t,x)};
            % @(t,x) (x(:,1) - x(:,2)) / (p.Tmax-t)};
            %@(t,x) (x(:,1) - x(:,2)) / t .* (t > 0)};

        % Second order term
        c.diffusion    = {@(t,x) 0.5*sigmaS^2*x(:,1).^2 .* (x(:,1) > p.xmin)  .*(x(:,1)<p.xmax);
            @(t,x) zeros(n(x),1);
            @(t,x) zeros(n(x),1);
            @(t,x) zeros(n(x),1)};

        c.g = @(t,x) max(x(:,2)-exp(-r*(p.Tmax-t))*K,0);
        c.f = @(t,x) zeros(n(x),1);

    end

    % Setup grid
    g = setupGrid(p);

end % function [p, c, g] = setupNonvariationalPDE(p, probLabel)
function ret = posConvection(t,x)
    ret = (x(:,1) - x(:,2)) / t;
    if t == 0
        ret = ones(size(x,1),1) .* ((x(:,1) > x(:,2)) - (x(:,1) < x(:,2)));
    end
end
