clear all
% Set bounds of domain
prob.xmin = 0.0;
prob.xmax = 1.0;
prob.ymin = 0.0;
prob.ymax = 1.0;

prob.domain = 'UnitSquare';
prob.dim = 2;
prob.DirichletBoundary = @(onBoundary,x) onBoundary;

% p.domain = 'Diamond';
[prob.geom,prob.gd] = getDomain(prob.domain, prob.xmin, prob.xmax, prob.ymin, prob.ymax);

% grid = setupUniformGrid(prob, 2);
% p = grid.p;
% e = grid.e;
% t = grid.t;

[p,e,t] = initmesh(prob.geom,'hmax',0.3);
grid.p = p;
grid.e = e;
grid.t = t;
grid.p_type = zeros(1,size(grid.p,2));
grid.p_type(unique(grid.e(1:2,:))) = 1;
grid.p_type(grid.e(1,grid.e(3,:)==0)) = 2;
grid.N = size(grid.p,2);

sfigure(1); clf; hold on
plotTriangulation(grid);
% return
corners = find(grid.p_type==2);
n_corners = length(corners);
idx_corners = 1:n_corners;

% Iniatialize V with corners
V = grid.p(:,corners)';

A = sort([grid.t([1,2],:), grid.t([2,3],:), grid.t([3,1],:)]);
[E,IA,IC] = unique(A','rows'); E = E';
n_edges = size(E,2);
idx_edges = n_corners + (1:n_edges);
% return
% Append midpoints of edges to vertex list
V = [V; 0.5*(grid.p(:,E(1,:))+grid.p(:,E(2,:)))'];

% Append centroids
centroids = 1/3*(grid.p(:,grid.t(1,:)) + grid.p(:,grid.t(2,:)) + grid.p(:,grid.t(3,:)));
n_centroids = size(centroids,2);
idx_centroids = n_corners + n_edges + (1:n_centroids);

V = [V;centroids'];
for i = 1:grid.N

    % Find containing triangles
    Ti = find(any(grid.t(1:3,:)==i));

    % Find containing edges
    Ei = find(any(E==i));

    if grid.p_type(i) == 2
        Ci = find(corners==i);
        v_idx = [idx_corners(Ci),idx_centroids(Ti),idx_edges(Ei)];
    else
        v_idx = [idx_centroids(Ti),idx_edges(Ei)];
    end
    vtx = V(v_idx,:);
    v_center = mean(vtx);
    v_0 = vtx - v_center;
    alpha = atan2(v_0(:,2),v_0(:,1));
    [~,order]=sort(alpha);
    C{i} = v_idx(order);
end
grid.V = V;
grid.C = C;
plotMesh(grid);
