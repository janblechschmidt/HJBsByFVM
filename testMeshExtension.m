p.dim = 2;
p.domain = 'UnitSquare';
p.mesh_type = 'regular';
% p.hmax = 0.5;
% p.hmax = 0.33;
p.hmax = 0.1;
p.geom = getDomain(p.domain);

% Set bounds of domain
p.xmin = 0.0;
p.xmax = 1.0;
p.ymin = 0.0;
p.ymax = 1.0;

p.control_dim = 1;

% Setup grid
g = setupGrid(p);
sfigure(1); clf, hold on
ax = plotTriangulation(g);
p.hmax = p.hmax/2;
gfine = setupGrid(p);
sfigure(2); clf, hold on
ax = plotTriangulation(gfine);

n = g.N;
nfine = gfine.N;

% Define numbers of nodes as they are in the mesh
I = reshape(1:n, sqrt(n), sqrt(n));
Ifine = reshape(1:nfine, sqrt(nfine), sqrt(nfine));
col = Ifine(1:2:end,1:2:end);
row = I;
c2f = sparse(row(:), col(:), ones(n,1), n, nfine);

% Identify the previous nodes with the new ones
% i = 1:n
% j = floor((i-1)/sqrt(n))
% k = mod(i-1,sqrt(n))
% idx_new = (k*2+1)+j*2*sqrt(nfine)
% c2f1 = sparse(i, idx_new, ones(n,1), n, nfine)
% keyboard

% Find the nodes that are interpolated by values above and below
col = Ifine(2:2:end,1:2:end); 
row1 = I(1:end-1,:);
row2 = I(2:end,:);
c2f = c2f + sparse([row1(:);row2(:)],[col(:);col(:)],0.5*ones(2*(n-sqrt(n)),1),n, nfine);

% Find the nodes that are interpolated by values left and right
col = Ifine(1:2:end,2:2:end);
row1 = I(:,1:end-1);
row2 = I(:,2:end);

c2f = c2f + sparse([row1(:);row2(:)],[col(:);col(:)],0.5*ones(2*(n-sqrt(n)),1),n, nfine);

% Find the nodes that are interpolated by diagonal values
col = Ifine(2:2:end,2:2:end);
row1 = I(1:end-1,1:end-1);
row2 = I(2:end,2:end);

c2f = c2f + sparse([row1(:);row2(:)],[col(:);col(:)],0.5*ones(2*((sqrt(n)-1)^2),1),n, nfine);

% The column sum should be one
max(abs(sum(c2f,1)-1))
abs(sum(sum(c2f)) - nfine)
% Now we define a function that we want to extend
f = @(x) sin(2*pi*x(:,1)) .* sin(2*pi*x(:,2));

sfigure(1); clf, hold on
fx = f(g.p');
plotOnTriangulation(g.p,g.t',fx)
fxfine = c2f'*fx;

sfigure(2); clf, hold on
plotOnTriangulation(gfine.p,gfine.t',fxfine)

