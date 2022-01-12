function u = control_MinimumArrivalTime(p, g, c, v, upwinding)
    helperFunctions;
    n = g.N;
    d = p.control_dim;
    EPS = 1e-8;
    if nargin < 5
        upwinding = 1;
    end
    if upwinding
        Aforw = g.D1forw;
        Aback = g.D1back;
    else
        Aforw = g.D1;
        Aback = g.D1;
    end

    u = zeros(n, d);
    for i = 1:d
        mu = @(u) c.convection{i}(0,0,u);
        Azero = @(u) spdiags(posPart(mu(u)),0,n,n) * Aforw{i} - ...
            spdiags(negPart(mu(u)),0,n,n) * Aback{i};
        val_func = @(u) Azero(u) * v + c.f(0,0,u);

        ustar_zero = zeros(n, d);
        ustar_forw = zeros(n, d);
        ustar_back = zeros(n, d);
        if p.alpha > 0
            ustar_forw(:, i) = (-Aforw{i} * v - p.beta) / p.alpha;
            ustar_back(:, i) = (-Aback{i} * v + p.beta) / p.alpha;
            % ustar_forw(:, i) = (-g.D1{i} * v - p.beta) / p.alpha;
            % ustar_back(:, i) = (-g.D1{i} * v + p.beta) / p.alpha;
        else
            ustar_forw(:, i) = p.Umax;
            ustar_back(:, i) = p.Umin;
        end
        % ustar_forw
        % keyboard
        % Truncate controls
        ustar_forw(ustar_forw < 0) = 0;
        ustar_back(ustar_back > 0) = 0;
        
        vforw = val_func(ustar_forw);
        vback = val_func(ustar_back);
        vzero = val_func(ustar_zero);
    
        V = [vzero, vforw, vback];
        % V(6,:)
        % [ustar_zero(6,i), ustar_forw(6,i), ustar_back(6,i)]
        [vmin, ind] = min(V, [], 2);

        % ind(abs(vmin - vforw)<EPS)=2;
        % ind(abs(vmin - vback)<EPS)=3;
        % ind(abs(vmin - vzero)<EPS)=1;
        u(:, i) = ustar_zero(:, i) .* (ind==1) + ustar_forw(:, i) .* (ind == 2) + ustar_back(:, i) .* (ind == 3);
    end
end
