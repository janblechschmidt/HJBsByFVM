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

grid = setupUniformGrid(prob, 8);
p = grid.p;
e = grid.e;
t = grid.t;
sfigure(1); clf; hold on
plotTriangulation(grid);

corners = find(grid.p_type==2);
n_corners = length(corners);
idx_corners = 1:n_corners;

% Iniatialize V with corners
V = grid.p(:,corners)';

A = sort([grid.t([1,2],:), grid.t([2,3],:), grid.t([3,1],:)]);
[E,IA,IC] = unique(A','rows'); E = E';
n_edges = size(E,2);
idx_edges = n_corners + (1:n_edges);

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
    Ti = find(any(grid.t==i));

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
    alpha = atan2(v_0(:,2),v_0(:,1));
    [~,order]=sort(alpha);
    C{i} = v_idx(order);
end
grid.V = V;
grid.C = C;
plotMesh(grid);
