% This file provides the coefficients for the energy storage problem
function c = setupCoefficients(p)
EPS = 1e-10;
n = @(x) size(x,1);

if strcmp( p.mode, 'MinimumArrivalTime_xcontrol' )

	c.DirichletBC = 1;
	c.alpha = p.alpha;   % Coefficient for L2 cost term
	c.beta = p.beta;    % Coefficient for L1 cost term
	valx = 0.5*p.sigmax^2;
	valy = 0.5*p.sigmay^2;
	valz = 0.5*p.sigmaz^2;
	valxz = 0.5 * p.corrxz * p.sigmax * p.sigmaz;
	valxy = 0.5 * p.corrxy * p.sigmax * p.sigmay;
	valyz = 0.5 * p.corryz * p.sigmay * p.sigmaz;
	
	if p.dim == 2
		c.diffusion    = {@(x) valx*ones(n(x),1);
		                  @(x) valxy*ones(n(x),1);
		                  @(x) valxy*ones(n(x),1);
		                  @(x) valy*ones(n(x),1)};
	
		c.convection   = {@(x,t,u) +u + p.mux *ones(n(x),1);
		                  @(x,t,u) p.muy*ones(n(x),1)};
	else
		if p.dim == 3
			c.diffusion    = {@(x) valx*ones(n(x),1);
			                  @(x) valxy*ones(n(x),1);
			                  @(x) valxz*ones(n(x),1);
			                  @(x) valxy*ones(n(x),1);
			                  @(x) valy*ones(n(x),1);
			                  @(x) valyz*ones(n(x),1);
			                  @(x) valxz*ones(n(x),1);
			                  @(x) valyz*ones(n(x),1);
			                  @(x) valz*ones(n(x),1)};
	
		c.convection   = {@(x,t,u) +u + p.mux *ones(n(x),1);
		                  @(x,t,u) p.muy*ones(n(x),1);
		                  @(x,t,u) p.muz*ones(n(x),1)};

		else
			error('Coefficients for dimension %i not implemented\n', p.dim)
		end % if p.dim == 3
	end % if p.dim == 2
	c.potential    = @(x) p.r * ones(n(x),1);
	% TODO: Changed sign of potential term
	c.finalTimeVal = @(x) zeros(n(x),1);
	c.boundaryVal  = @(x) zeros(n(x),1);
	c.f            = @(x,u) (1 + 0.5 * c.alpha * u.^2 + c.beta * abs(u));

end % if strcmp( p.mode, 'MinimumArrivalTime_xcontrol' )

if strcmp( p.mode, 'MinimumArrivalTime' )

	c.DirichletBC = 1;
	c.alpha = p.alpha;   % Coefficient for L2 cost term
	c.beta = p.beta;    % Coefficient for L1 cost term
	valx = 0.5*p.sigmax^2;
	valy = 0.5*p.sigmay^2;
	valxy = 0.5 * p.corrxy * p.sigmax * p.sigmay;
	
	if p.dim == 2
		c.diffusion    = {@(x) valx*ones(n(x),1);
											@(x) valxy*ones(n(x),1);
											@(x) valxy*ones(n(x),1);
											@(x) valy*ones(n(x),1)};

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
	
		c.convection   = {@(x,t,u) u(:,1) + p.mux *ones(n(x),1);
		                  @(x,t,u) u(:,2) + p.muy*ones(n(x),1)};
	end % if p.dim == 2

	c.potential    = @(x) p.r * ones(n(x),1);
	% TODO: Changed sign of potential term
	c.finalTimeVal = @(x) zeros(n(x),1);
	c.boundaryVal  = @(x) zeros(n(x),1);
	c.f            = @(x,u) (1 + 0.5 * c.alpha * sqrt(u(:,1).^2 + u(:,2).^2) + c.beta * (abs(u(:,1)) + abs(u(:,2))));

end % if strcmp( p.mode, 'MinimumArrivalTime' )

end % function setupCoefficients(p)
