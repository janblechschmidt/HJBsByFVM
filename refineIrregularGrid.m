function grid = refineIrregularGrid(prob, grid, it)
% This function has to update all grid-relevant quantities, both for the primal and dual mesh

if nargin < 3
    [p,e,t] = refinemesh(prob.geom, grid.p, grid.e, grid.t);
else
    [p,e,t] = refinemesh(prob.geom, grid.p, grid.e, grid.t, it);
end

grid.p = p;
grid.e = e;
grid.t = t;
% e = grid.e;
% p = grid.p;
% t = grid.t(1:3,:);
% np = size(grid.p,2);
% % x coordinates of new points along boundary
% px = 0.5*(p(1,e(1,:))+p(1,e(2,:)));
% py = 0.5*(p(2,e(1,:))+p(2,e(2,:)));
% new_np = size(e,2);
% new_boundary_points = [px;py];
% % new.p = [grid.p,[px;py]];

% % Compute new edges (they are not directed, but always have the lower index in the first row)
% % new.e = reshape([grid.e(1,:); (1:new_np)+np; grid.e(2,:); (1:new_np)+np],2,[]);
% t_new = [];
% keyboard
% for i = 1:grid.N_tri
%     ti = t(:,i);

%     e_ti_1 = sort(ti([1,2]));
%     e_ti_2 = sort(ti([2,3]));
%     e_ti_3 = sort(ti([3,1]));

%     idx_1 = find(all(grid.e==e_ti_1))+np;
%     idx_2 = find(all(grid.e==e_ti_2))+np;
%     idx_3 = find(all(grid.e==e_ti_3))+np;
    
%     % Append new triangles with the same orientation as the old ones
%     t_new = [t_new, [ti(1);idx_1;idx_3], [ti(2);idx_2;idx_1], [ti(3);idx_3;idx_2], [idx_1;idx_2;idx_3] ];
    
% end
% % new.t = t_new

% % volT = grid.muT*[.25,.25,.25,.25];
% % new.muT = reshape(volT',[],1);

% % keyboard

% % I need
% grid.p = [grid.p,[px;py]];
% grid.t = t_new;
% grid.e = reshape([grid.e(1,:); (1:new_np)+np; grid.e(2,:); (1:new_np)+np],2,[]);

if ~isfield(grid,'p_type') | length(grid.p_type)~=size(grid.p,2)
    grid.p_type = zeros(1,size(grid.p,2));
    grid.p_type(unique(grid.e(1:2,:))) = 1;
    grid.p_type(grid.e(1,grid.e(3,:)==0)) = 2;
end


grid.N = size(grid.p,2);          % Number of vertices == Number of control volumes
grid.N_tri = size(grid.t,2); % Number of triangles

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Compute dual grid  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch prob.dual_mesh
case 'voronoi'

    corners = find(grid.p_type==2);
    n_corners = length(corners);
    idx_corners = 1:n_corners;
    
    % Iniatialize V with corners
    V = grid.p(:,corners)';
    
    A = sort([grid.t([1,2],:), grid.t([2,3],:), grid.t([3,1],:)]);
    [E,IA,IC] = unique(A','rows');
    E = E';
    n_edges = size(E,2);
    idx_edges = n_corners + (1:n_edges);
    grid.E = E;
    
    % Append midpoints of edges to vertex list
    v_edge = 0.5*(grid.p(:,E(1,:))+grid.p(:,E(2,:)));
    % V = [V; v_edge'];
    
    % Append centers
    e_a = grid.t([1,2],:);
    e_b = grid.t([2,3],:);
    a = .5*(grid.p(:,e_a(1,:)) + grid.p(:,e_a(2,:)));
    b = .5*(grid.p(:,e_b(1,:)) + grid.p(:,e_b(2,:)));
    a_vec = grid.p(:,e_a(1,:)) - grid.p(:,e_a(2,:));
    b_vec = grid.p(:,e_b(1,:)) - grid.p(:,e_b(2,:));
    a_orth = [a_vec(2,:);-a_vec(1,:)];
    b_orth = [b_vec(2,:);-b_vec(1,:)];
    det_V = -a_orth(1,:).*b_orth(2,:)+b_orth(1,:).*a_orth(2,:);
    lambda = 1./det_V.*( -b_orth(2,:).*(b(1,:)-a(1,:)) + b_orth(1,:).*(b(2,:)-a(2,:)) );
    v_mid = a + lambda.*a_orth;
    
    % It might occur that the new vertices are already added as edge midpoints
    [v_new, ia, ic] =  unique([v_edge,v_mid]','rows','stable');
    
    idx_centers = n_corners + ic(n_edges+1:end)';
    
    V = [V; v_new];
    
    for i = 1:grid.N
    
        % Find containing triangles
        Ti = find(any(grid.t(1:3,:)==i));
    
        % Find containing edges
        Ei = find(any(E==i));
    
        if grid.p_type(i) == 2
            Ci = find(corners==i);
            v_idx = unique([idx_corners(Ci),idx_centers(Ti),idx_edges(Ei)]);
        else
            v_idx = unique([idx_centers(Ti),idx_edges(Ei)]);
        end
        vtx = V(v_idx,:);
        v_center = mean(vtx);
        v_0 = vtx - v_center;
        % Ensure clockwise ordering
        alpha = atan2(v_0(:,2),v_0(:,1));
        [~,order]=sort(alpha);
        C{i} = v_idx(order);

    end
    grid.C = C;
    grid.V = V;

    % grid.dt = delaunayTriangulation;
    % grid.dt.Points = p';
    % % This case is currently treated using the MATLAB routine
    % [V,C] = voronoiDiagram(grid.dt);
    % 
    % % DEBUG: Plot unbounded Voronoi cells
    % % voronoi(grid.dt)
    % 
    % % Out of some reason, voronoiDiagram is not that much stable,
    % % in particular in doesn't recognize unbounded cells reliably
    % inf_like = find(sum(V.^2,2)> (10000+(prob.xmax-prob.xmin).*(prob.ymax-prob.ymin)));
    % if length(inf_like) > 1
    %     % Delete multiple infs
    %     for i=1:grid.N
    %         % if i==5
    %         %     keyboard
    %         % end
    % 
    %         [int,ia,ib] = intersect(inf_like,C{i});
    %         % if length(int)>2
    %         %     warning('More than three infs in voronoi diagramm\n')
    %         %     keyboard
    %         % end
    %         if length(int) > 0
    %             if int(1)==1 & length(int) > 1 % first entry is (int,int)
    %                 C{i}(ib(ib~=1))=[];
    %                 % C{i}(C{i}==int(2))=[];
    %             else
    %                 if length(int) == 1
    %                     j = find(int==C{i});
    %                     C{i} = [1,C{i}((j+1):end),C{i}(1:(j-1))];
    %                 else
    %                     % if i==8
    %                     %     keyboard
    %                     % end
    %                     % Find occurence of first inf_like vertex
    %                     j = min(ib);
    %                     % j = find(int(1)==C{i});
    %                     % Delete the other ones
    %                     C{i}(ib(ib~=j)) = [];
    %                     C{i} = [1,C{i}((j+1):end),C{i}(1:(j-1))];
    %                 end
    %             end
    %         end
    %     end
    %    
    %     % Now, all the indices have changed, let's correct this
    %     M=size(V,1);
    %     correction = (sum((1:M)'>=repmat(inf_like(2:end)',M,1),2));
    %     for i=1:grid.N
    %         C{i} = C{i} - correction(C{i})';
    %     end
    %     
    %     % Now, we can safely remove the unnecessary inf duplicates
    %     V(inf_like(2:end),:) = [];
    % end
    % 
    % % initMesh sometimes generates meshes which contain Voronoi vertices
    % % outside of the domain. The following function finds these points
    % % and corrects the primal mesh in order to correct the dual mesh.
    % DV = computeDistanceToBoundary(prob.geom, V(2:end,:)');
    % Vidx = find(any(DV>1e-10,2));
    % 
    % % If there are Voronoi vertices outside of the domain, we have to add 
    % % the projections of these points onto the domain as a primal node
    % 
    % if ~isempty(Vidx)
    % 	warning('Triangulation will be adapted to be suited for perpendicular bisector finite volumes');
    %     warning('TODO: determineClosestBndValue has to be adapted for usage on non-rectangular domains\n')
    %     keyboard
    % 	 % Find Voronoi vertices outside of the domain and add them as points to p
    % 	 for v=Vidx'
    % 		[v_at_bnd] = determineClosestBndValue(prob, V(v+1,:));
    % 		p(:,end+1) = v_at_bnd{1}';
    % 	 end
    % 	grid.N = size(p,2); % Number of vertices == Number of control volumes
    % 	grid.dt = delaunayTriangulation(p(1,:)',p(2,:)');
    % 	grid.N_tri = size(grid.dt.ConnectivityList,1);
    %     grid.t = [grid.dt.ConnectivityList';ones(1,grid.N_tri)];
    % 	grid.p = p;
    % 
    % 	% Setup dual grid
    % 
    % 	[V,C] = voronoiDiagram(grid.dt);
    %     DV = computeDistanceToBoundary(prob.geom, V(2:end,:)');
    % end
    % 
    % % The function voronoiDiagram gives unbounded Voronoi Cells, but we need bounded ones.
    % [V, C, inside] = computeBoundedVoronoiCells(prob, grid, V, C, DV);

case 'centroid'
    % keyboard

    grid = setupDualCentroidGrid(grid, prob);
end
inside = (grid.p_type == 0)';
grid.DomainBoundary = find(~inside);
grid.DirichletBoundary = find(prob.DirichletBoundary(~inside, p'));
grid.OutflowBoundary = setdiff(grid.DomainBoundary,grid.DirichletBoundary);
grid.inner = setdiff(1:grid.N, grid.DirichletBoundary);

% grid.E = E;
% grid.N_inner = sum(inside);
% grid.C = C;
grid = computeDualMeshQuantities(grid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Compute mesh specific quantities  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grid = setupGridBoundaries(grid, prob);
grid = computeMeshQuantities(grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Precompute matrices for differential operators  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

grid = setupZerothOrderMatrix(grid);
grid = setupFirstOrderMatrix(grid);
grid = setupSecondOrderMatrix(grid);
end % function refineIrregularGrid(prob, grid)
