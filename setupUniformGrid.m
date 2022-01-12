function grid = setupUniformGrid(prob, mesh_N)


    % Problem dimension is also relevant for grid structure
    grid.dim = prob.dim;
    grid.mesh_N = mesh_N;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute primal grid %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ~isfield(prob,'grid_type')
        grid.type = 'left';
    else
        grid.type = prob.grid_type;
    end

    % Determine vertices in primal mesh
    switch prob.dim
        case 1

            x = linspace(prob.xmin, prob.xmax, mesh_N + 1);

            % List of vertex coordinates
            grid.p = x;

            % List of cells
            grid.t = [1:mesh_N;2:mesh_N+1];


        case 2
            x=linspace(prob.xmin, prob.xmax, mesh_N + 1);
            y=linspace(prob.ymin, prob.ymax, mesh_N + 1);
            [X,Y] = meshgrid(x,y);

            grid.p = [X(:)'; Y(:)'];
            hx = x(2)-x(1);
            hy = y(2)-y(1);

            tmp=zeros(mesh_N+1,1);
            tmp(1) = 1; tmp(end) = 1;
            tmp = tmp+tmp';
            grid.p_type = tmp(:);

            switch grid.type
                case 'left'

                    % Set list of triangles
                    T = reshape(1:(mesh_N+1)^2,mesh_N+1,mesh_N+1)';
                    T1 = reshape(T(1:end-1,1:end-1),[],1);
                    T2 = reshape(T(1:end-1,2:end  ),[],1);
                    T3 = reshape(T(2:end  ,2:end  ),[],1);
                    U1 = reshape(T(1:end-1,1:end-1),[],1);
                    U2 = reshape(T(2:end  ,2:end  ),[],1);
                    U3 = reshape(T(2:end  ,1:end-1),[],1);
                    grid.t = [T1, T2, T3; U1, U2, U3]';

                case 'right'

                    % Set list of triangles
                    T = reshape(1:(mesh_N+1)^2,mesh_N+1,mesh_N+1)';
                    % Flip T to make triangle assembling easier
                    T = fliplr(T);
                    T1 = reshape(T(1:end-1,1:end-1),[],1);
                    T2 = reshape(T(2:end  ,1:end-1),[],1);
                    T3 = reshape(T(2:end  ,2:end  ),[],1);
                    U1 = reshape(T(1:end-1,1:end-1),[],1);
                    U2 = reshape(T(2:end  ,2:end  ),[],1);
                    U3 = reshape(T(1:end-1,2:end  ),[],1);
                    grid.t = [T1, T2, T3; U1, U2, U3]';
                    % Flip T back to original state
                    T = fliplr(T);

            end
            % Set list of edges
            tmp = linspace(0, 1, mesh_N + 1);
            e1 = [T(1:end-1,1)';
                T(2:end,1)';
                tmp(1:end-1);
                tmp(2:end);
                2*ones(1,mesh_N);
                ones(1,mesh_N);
                zeros(1,mesh_N)];
            e2 = [T(end,1:end-1);
                T(end,2:end);
                tmp(1:end-1);
                tmp(2:end);
                3*ones(1,mesh_N);
                ones(1,mesh_N);
                zeros(1,mesh_N)];
            e3 = [fliplr(T(2:end,end)');
                fliplr(T(1:end-1,end)');
                tmp(1:end-1);
                tmp(2:end);
                4*ones(1,mesh_N);
                ones(1,mesh_N);
                zeros(1,mesh_N)];
            e4 = [fliplr(T(1,2:end));
                fliplr(T(1,1:end-1));
                tmp(1:end-1);
                tmp(2:end);
                1*ones(1,mesh_N);
                ones(1,mesh_N);
                zeros(1,mesh_N)];
            grid.e = [e4,e1,e2,e3];

    end

    
    grid = setupGridBoundaries(grid, prob);

    % Assign primal quantities to grid structure
    grid.N = size(grid.p,2);          % Number of vertices == Number of control volumes
    grid.N_tri = size(grid.t,2); % Number of triangles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Compute dual grid  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch prob.dim
        case 1
            hx = grid.p(2) - grid.p(1);
            grid.muGamma = spdiags(ones(grid.N,2),[-1,1],grid.N,grid.N);
            grid.Ni = spdiags(ones(grid.N,2),[-1,1],grid.N,grid.N);
            grid.Nij{1} = spdiags([-ones(grid.N,1),ones(grid.N,1)],[-1,1],grid.N,grid.N);
            grid.inside = ones(grid.N,1);
            grid.inside(1) = 0;
            grid.inside(grid.N) = 0;
            V = [grid.p(1), 0.5*(grid.p(2:end) + grid.p(1:end-1)), grid.p(end)];

            muOmega = hx*ones(grid.N,1);
            muOmega(1) = hx/2;
            muOmega(end) = hx/2;

            grid.muOmega = muOmega;

            for i = 1:mesh_N+1
                C{i} = [i,i+1];
            end % for i = 1:mesh_N

            grid.orientationT = ones(grid.N_tri,1);
            grid.muT = hx*ones(grid.N_tri,1);

            grid.C = C;
        case 2
            grid.orientationT = -ones(grid.N_tri,1);
            if strcmpi(prob.domain, 'Diamond')
                scale = 1/sqrt(2);
            else
                scale = 1;
            end
            grid.muT = scale^2*0.5*hx*hy*ones(grid.N_tri,1);
            switch prob.dual_mesh
                case 'voronoi'
                    % The regular case is much more easier, and needs no real computation
                    xmid = 0.5*(x(2:end) + x(1:end-1));
                    ymid = 0.5*(y(2:end) + y(1:end-1));
        
                    xmid = [prob.xmin, xmid, prob.xmax];
                    ymid = [prob.ymin, ymid, prob.ymax];
        
                    [Xmid, Ymid] = meshgrid(xmid, ymid);
                    V = [Xmid(:), Ymid(:)];
        
                    if strcmpi(prob.domain, 'Diamond')
                        % Set up rotation matrix (with scaling)
                        alpha = 45/180*pi;
                        PREC = 14;
                        scale = 1/sqrt(2);
                        S = [0.5*(prob.xmin+prob.xmax);0.5*(prob.ymin+prob.ymax)];
                        R = scale*[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
                        % Rotate the vertices of the dual mesh
                        grid.V = round((V-S')*R'+S',14);
                        % as well as the nodes of the primal mesh
                        grid.p = round(R*(grid.p-S)+S,14);
                    else
                        scale = 1;
                        grid.V = V;
                    end
        
                    % There are nine different cell categories to distinguish
                    globidx = flipud(reshape(1:(mesh_N+1)^2,mesh_N+1,mesh_N+1));
                    % - lower left corner
                    llc = reshape(globidx(end,1),[],1);
                    % - left boundary
                    leb = reshape(globidx(2:end-1,1),[],1);
                    % - upper left corner
                    ulc = reshape(globidx(1,1),[],1);
                    % - lower boundary
                    lob = reshape(globidx(end,2:end-1),[],1);
                    % - inner cells
                    inc = reshape(globidx(2:end-1,2:end-1),[],1);
                    % - upper boundary
                    upb = reshape(globidx(1,2:end-1),[],1);
                    % - lower right corner
                    lrc = reshape(globidx(end,end),[],1);
                    % - right boundary
                    rib = reshape(globidx(2:end-1,end),[],1);
                    % - upper right boundary
                    urc = reshape(globidx(1,end),[],1);
        
                    % Set volume of Voronoi cells
                    muOmega = zeros(grid.N,1);
        
                    muOmega_inner  = hx * hy;
                    muOmega_bnd    = muOmega_inner / 2;
                    muOmega_corner = muOmega_inner / 4;
        
                    muOmega(llc) = muOmega_corner;
                    muOmega(ulc) = muOmega_corner;
                    muOmega(lrc) = muOmega_corner;
                    muOmega(urc) = muOmega_corner;
        
                    muOmega(leb) = muOmega_bnd;
                    muOmega(lob) = muOmega_bnd;
                    muOmega(rib) = muOmega_bnd;
                    muOmega(upb) = muOmega_bnd;
        
                    muOmega(inc) = muOmega_inner;
        
                    grid.muOmega = scale^2*muOmega;
        
                    % Set inside vector
                    grid.inside = zeros(grid.N,1);
        
                    grid.inside(inc) = 1;
        
                    % Set lengths of edges between Voronoi cells
                    L = [llc, llc+1, 0.5*hx;
                        llc, llc+mesh_N+1, 0.5*hy;
                        leb, leb-1, 0.5*hx*ones(mesh_N-1,1);
                        leb, leb+mesh_N+1, hy*ones(mesh_N-1,1);
                        leb, leb+1, 0.5*hx*ones(mesh_N-1,1);
                        ulc, ulc-1, 0.5*hx;
                        ulc, ulc+mesh_N+1, 0.5*hy;
                        lob, lob-mesh_N-1, 0.5*hy*ones(mesh_N-1,1);
                        lob, lob+1, hx*ones(mesh_N-1,1);
                        lob, lob+mesh_N+1, 0.5*hy*ones(mesh_N-1,1);
                        inc, inc-mesh_N-1, hy*ones((mesh_N-1)^2,1);
                        inc, inc-1, hx*ones((mesh_N-1)^2,1);
                        inc, inc+1, hx*ones((mesh_N-1)^2,1);
                        inc, inc+mesh_N+1, hy*ones((mesh_N-1)^2,1);
                        upb, upb-mesh_N-1, 0.5*hy*ones(mesh_N-1,1);
                        upb, upb-1, hx*ones(mesh_N-1,1);
                        upb, upb+mesh_N+1, 0.5*hy*ones(mesh_N-1,1);
                        lrc, lrc-mesh_N-1, 0.5*hy;
                        lrc, lrc+1, 0.5*hx;
                        rib, rib-mesh_N-1, hy*ones(mesh_N-1,1);
                        rib, rib-1, 0.5*hx*ones(mesh_N-1,1);
                        rib, rib+1, 0.5*hx*ones(mesh_N-1,1);
                        urc, urc-mesh_N-1, 0.5*hy;
                        urc, urc-1, 0.5*hx];
        
                    grid.muGamma = scale*sparse(L(:,1), L(:,2), L(:,3),grid.N,grid.N);
        
                    % Set up matrix that stores information that was previously in Ni
                    grid.Ni = sparse(L(:,1), L(:,2), ones(size(L(:,3))),grid.N,grid.N);
        
                    nij = grid.p(:,L(:,2))-grid.p(:,L(:,1));
                    muij = sqrt(sum(nij.^2,1));
        
                    grid.Nij{1} = sparse(L(:,1), L(:,2), nij(1,:)./muij, grid.N, grid.N);
                    grid.Nij{2} = sparse(L(:,1), L(:,2), nij(2,:)./muij, grid.N, grid.N);
                    grid.Nuij{1} = grid.muGamma.*grid.Nij{1};
                    grid.Nuij{2} = grid.muGamma.*grid.Nij{2};
        
                    for i = 1:(mesh_N+1)
                        for j = 1:(mesh_N+1)
        
                            % Determine global index
                            idx = (i-1)*(mesh_N+1)+j;
        
                            % Determine the indices of the four corners of the Voronoi cell with global index idx
                            C{idx} = [(i-1)*(mesh_N+2)+j, i*(mesh_N+2)+j, i*(mesh_N+2)+1+j, (i-1)*(mesh_N+2)+j+1];
        
                        end % for j = 1:mesh_N
                    end % for i = 1:mesh_N
                    grid.C = C;

                case 'centroid'
                    grid = setupDualCentroidGrid(grid);
                    grid = computeDualMeshQuantities(grid);
            end
    end

    grid.t(4,:) = 1;

    % Assign some dual quantities to grid structure


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Compute mesh specific quantities  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine and boundary vertices
    grid.DomainBoundary = find(grid.lowbnd{1} | grid.lowbnd{2} | grid.upbnd{1} | grid.upbnd{2});
    % grid.bnd = find(prob.DirichletBoundary((grid.lowbnd{1} | grid.lowbnd{2} | grid.upbnd{1} | grid.upbnd{2})', grid.p'));
    grid.DirichletBoundary = grid.DomainBoundary(find(prob.DirichletBoundary(grid.DomainBoundary', grid.p(:,grid.DomainBoundary)')));
    grid.OutflowBoundary = setdiff(grid.DomainBoundary, grid.DirichletBoundary);
    grid.InnerIndices = setdiff(1:(mesh_N+1)^prob.dim, grid.DomainBoundary);
    grid.inner = union(grid.InnerIndices,grid.OutflowBoundary);

    % grid.bnd = find(prob.DirichletBoundary(~grid.inside, grid.p'));
    % grid.inner = setdiff(1:(mesh_N+1)^prob.dim, grid.bnd);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Precompute matrices for differential operators  %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    grid = setupZerothOrderMatrix(grid);
    grid = setupFirstOrderMatrix(grid);
    grid = setupSecondOrderMatrix(grid);

    % grid2 = computeMeshQuantities(grid);
    % max(max(abs(grid.muT-grid2.muT)))
    % max(max(abs(grid.muOmega-grid2.muOmega)))
    % max(max(abs(grid.muGamma-grid2.muGamma)))
    % max(max(abs(grid2.Ny.*grid.Ni-grid.Ny)))
