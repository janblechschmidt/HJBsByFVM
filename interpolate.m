function [vx] = interpolate(X, v, grid)
X = X';
lowbnds = min(X,[],2);
upbnds = max(X,[],2);

[m,n] = size(X);
X = reshape(X,m,1,n);
fprintf('Interpolate grid...\n')

x0 = grid.p(:,grid.t(1,:));
x1 = grid.p(:,grid.t(2,:));
x2 = grid.p(:,grid.t(3,:));
EPS=1e-10;

% Determine relevant triangles to improve performance
Threshold = 1e5;
rad = 8*full(max(max(grid.muGamma)));
if grid.N_tri * n > Threshold
    xm = 1./3*(x0+x1+x2);
    tri_rel = find(all(xm > lowbnds-rad & xm < upbnds+rad));
else
    tri_rel = 1:grid.N_tri;
end

v1 = x1(:,tri_rel)-x0(:,tri_rel);
v2 = x2(:,tri_rel)-x0(:,tri_rel);
v0 = X-x0(:,tri_rel);

% The following could be precomputed and stored during the computation of triangle size
    % detVinv = 1./(v1(1,:).*v2(2,:)-v2(1,:).*v1(2,:));
detVinv = (1./(grid.muT(tri_rel)*2).*grid.orientationT(tri_rel))';
% If n and grid.N_tri is large, we have to split the computation
if length(tri_rel) * n > Threshold
    ns = floor(length(tri_rel) * n / Threshold);
    ns = max(round(n/ns),1);
else
    % ns = 1;
    ns = n;
end
startidx=1;
vx = zeros(1,n);
done = 0;
while ~done
    fprintf('.')
    if startidx + ns < n
        idx = startidx:(startidx+ns);
        startidx = startidx+ns+1;
    else
        idx = startidx:n;
        done = 1;
    end
    nidx = length(idx);
    lambda = detVinv.*(v2(2,:).*v0(1,:,idx)-v2(1,:).*v0(2,:,idx));
    mu = detVinv.*(-v1(2,:).*v0(1,:,idx)+v1(1,:).*v0(2,:,idx));
    
    % Find the elements containing the points
    tmp =  lambda>=0-EPS & lambda<=1+EPS & mu>=0-EPS & mu<=1+EPS & lambda+mu<=1+EPS;
    tmp = reshape(tmp,[],nidx)';
    i = 1:nidx;
    % size(tmp)
    % tmp
    % idx
    % X(:,idx)
    j = arrayfun(@(x) find(tmp(x,:),1),i);
    
    % [~,j] = max(tmp,[],2);
    % j = j';
    globidx = sub2ind(size(lambda),ones(size(i)),j,i);
    tx = grid.t(:,tri_rel(j));
    s3 = mu(globidx);
    s2 = lambda(globidx);
    s1 = (1-s2-s3);
    
    vx(idx) = s1.*v(tx(1,:))'+s2.*v(tx(2,:))'+s3.*v(tx(3,:))';
    % max(abs(sum([1-mu(globidx)'-lambda(globidx)', lambda(globidx)', mu(globidx)'],2)-1))
end
fprintf('\n')
end % function vx = interpolate(X, grid)
