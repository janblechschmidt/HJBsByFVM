function plotOnTriangulation(p, c, v)
% INPUT: p - matrix of size 2 x N of coordinates in R^2
%        c - connectivity of triangles of size M x 3
%        v - vector of size N (CG1) or M (DG0) of values of function
M = size(c,1); % Number of triangles
[dim,N] = size(p); % Number of vertices

switch dim
case 1
    plot(p,v,'r+-')
case 2
    % Prepare x and y coordinates
    X = zeros(M,4);
    Y = zeros(M,4);
    Z = zeros(M,4);
    
    if size(v,1) == N
    	for j = 1:M
    		tri = c(j,1:3);
    		idx = [tri, tri(1)];
    		X(j,:) = p(1,idx);
    		Y(j,:) = p(2,idx);
    		Z(j,:) = v(idx);
    	end
    end
    
    if size(v,1) == M
    	for j = 1:M
    		tri = c(j,:);
    		idx = [tri, tri(1)];
    		X(j,:) = p(1,idx);
    		Y(j,:) = p(2,idx);
    		Z(j,:) = v(j);
    	end
    end
    
    % Plot cells
    % azsurf = fill3(X', Y', Z', Z','EdgeAlpha',0.2,'LineWidth',0.01);
    azsurf = fill3(X', Y', Z', Z','EdgeColor','none','LineWidth',0.01);
    colormap viridis;
    alpha(azsurf,0.9);
end % switch dim

end % function plot_f_on_tri(p,c,v)
