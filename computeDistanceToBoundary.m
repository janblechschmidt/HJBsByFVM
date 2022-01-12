function d = computeDistanceToBoundary(geom, X)
m = size(X,2);
n = size(geom,2);
d = zeros(m,n);

for ei = 1:n
    p1 = geom([2,4], ei);
    p2 = geom([3,5], ei);
    u = p2-p1;
    uN = u / sqrt(u'*u);
    uON = [uN(2);-uN(1)];
    lambda = uN' * (X-p1);
    projx = p1 + lambda .* uN;
    d(:,ei) = sum((X-p1).*uON)';
    % distx = sum((X-p1).*(X-projx),1);
    % dX = sum(X.*(X-projx),1)-distx;
    % d(:,ei) = dX';
end

end
