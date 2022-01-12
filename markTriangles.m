function markedTriangles = markTriangles(v, grid)

markedTriangles = [];
markedPoints = [];

% Compute edge-based error estimator for each triangle
eta = zeros(grid.N_tri,3);

% Compute some triangle related values

% First dimension gives triangle index
grad = zeros(grid.N_tri,2);
el   = zeros(grid.N_tri,2,3);
nv   = zeros(grid.N_tri,2,3);
le   = zeros(grid.N_tri,3);

for i = 1:grid.N_tri

    tri = grid.t(:,i)';
	% tri = grid.cl(i,:);
	coords = grid.p(:,tri);
	v_local = v(tri);

	% Get midpoints of cell edges
	el(i,:,1) = mean(coords(:,[1,2]),2);
	el(i,:,2) = mean(coords(:,[2,3]),2);
	el(i,:,3) = mean(coords(:,[3,1]),2);

	% Compute length of edges
	le(i,1) = sqrt(sum(diff(coords(:,[1,2]),[],2).^2));
	le(i,2) = sqrt(sum(diff(coords(:,[2,3]),[],2).^2));
	le(i,3) = sqrt(sum(diff(coords(:,[3,1]),[],2).^2));

	% Get normal vector on cell edges
	nv(i,:,1) = diff(coords(:,[1,2]),[],2) / le(i,1);
	nv(i,:,2) = diff(coords(:,[2,3]),[],2) / le(i,2);
	nv(i,:,3) = diff(coords(:,[3,1]),[],2) / le(i,3);
	tmp = nv(i,1,:);
	nv(i,1,:) = -nv(i,2,:);
	nv(i,2,:) = tmp;
	if nv(i,:,1)*nv(i,:,2)' < 0
		nv(i,:,:) = -nv(i,:,:);
	end

	grad(i,:) = gradTri(coords,v_local);

end % for i = 1:grid.N_tri

% Plot the gradient on each triangle
sfigure(4); clf; hold on
subplot(1,2,1);% clf; hold on
plotOnTriangulation(grid.p, grid.t', grad(:,1))
subplot(1,2,2);% clf; hold on
plotOnTriangulation(grid.p, grid.t', grad(:,2))
% plotOnTriangleGrid(grid, grad, ['v_x'; 'v_y'])


edge_list = [1,2 ; 2,3 ; 3,1];
for i = 1:grid.N_tri

    tri = grid.t(:,i)';
	% tri = grid.cl(i,:);

	grad_tri  = grad(i,:)';

    keyboard
    b > find(grid.Ni(i,:)>0)
	% Find indices of neighbors of triangle i
	b = [any(grid.cl==tri(1),2), any(grid.cl==tri(2),2), any(grid.cl==tri(3),2)];
	b = [any(grid.t==tri(1),1), any(grid.t==tri(2),1), any(grid.t==tri(3),1)];

	for j = 1:3
		edge = edge_list(j,:);
		n_idx =  setdiff(find(b(:,edge(1)) & b(:,edge(2))),i);
		if ~isempty(n_idx)
			grad_n    = grad(n_idx,:)';
			keyboard
			grad_diff = grad_n - grad_tri;
			nv_tri    = nv(i,:,j)';
			eta(i,j) = le(i,j)^0.5 * abs( grad_diff' * nv_tri );
		else
			eta(i,j) = 0;
		end
	end
% 	n_idx(1) = setdiff(find(b(:,1) & b(:,2)),i);
% 	n_idx(2) = setdiff(find(b(:,2) & b(:,3)),i);
% 	n_idx(3) = setdiff(find(b(:,3) & b(:,1)),i);

end

eta_ell = max(eta,[],2);
eta_max = max(eta_ell);

[eta_sort,rearrange_idx] = sort(eta_ell,'descend');
% markedTriangles = rearrange_idx(cumsum(eta_sort) < opt.refinement_ratio *sum(eta_ell));

fprintf('||eta||_inf = %8.4f\n', eta_max);

% Plot maximum error estimate on triangles
sfigure(3); clf; hold on
plotOnTriangleGrid(grid, eta_ell)
xlabel('x'), ylabel('y')
title('Cellwise error estimate');
view([0,90])

% Plot point-wise error estimate
sfigure(5); clf; hold on
plotTriangulation(grid.dt);
X = el(:,1,:);
Y = el(:,2,:);
Z = eta;
C = Z / max(max(Z));
scatter3(X(:),Y(:),Z(:),30*ones(size(Z(:))),C(:)*[1 0 0],'o','filled');
title('Pointwise error estimate');
% view([50,20])
view([0,90])

% Alternative approach to give coordinates of new points instead of cells
[eta_sort,rearrange_idx] = sort(eta(:),'descend');
eta_ell = rearrange_idx(cumsum(eta_sort) < opt.refinement_ratio *sum(eta(:)));
[I,J,K]=ind2sub([size(el,1),1,size(el,3)],eta_ell);
ex = el(sub2ind([size(el)], I, ones(size(I)), K));
ey = el(sub2ind([size(el)], I, 2*ones(size(I)), K));
markedPoints = unique([ex,ey],'rows');
pause
num_ref = num_ref + 1
if num_ref > opt.number_of_refinements
	ref_done = 1;
end


% New strategy to evaluate the error estimator
% First step: extract the edges from triangle data structure
% 
% for i = 1:grid.N_tri
% 	tri = grid.cl(i,:);
% 	coords = grid.p(:,tri);
% 	v_local = v(tri);
% 	xij = zeros(2,3);
% 	muij = zeros(1,3);
% 	xij(1:2,1) = mean(coords(:,[1,2]),2);
% 	xij(1:2,2) = mean(coords(:,[2,3]),2);
% 	xij(1:2,3) = mean(coords(:,[3,1]),2);
% 
% 	% Compute outer-normal vectors and lengths of edges
% 	nv(:,1) = diff(coords(:,[1,2]),[],2);
% 	nv(:,2) = diff(coords(:,[2,3]),[],2);
% 	nv(:,3) = diff(coords(:,[3,1]),[],2);
% 
% 	muij(1) = sqrt(sum(nv(:,1).^2));
% 	muij(2) = sqrt(sum(nv(:,2).^2));
% 	muij(3) = sqrt(sum(nv(:,3).^2));
% 
% 	nv(:,1) = nv(:,1) / muij(1);
% 	nv(:,2) = nv(:,2) / muij(2);
% 	nv(:,3) = nv(:,3) / muij(3);
% 
% 	nv(1,:) = -nv(1,:);
% 	if nv(:,1) * nv(:,2)' > 0
% 		nv = -nv;
% 	end
% 	keyboard
% 
% end % for i = 1:grid.N_tri




