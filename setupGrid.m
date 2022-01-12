function [grid] = setupGrid(prob, pnts, e, t)

fprintf('Setup grid...\n')
if nargin < 2
    switch prob.mesh_type
    case 'regular'
    
        % Determine number of equidistant intervals to get
        % a maximum mesh size of approximately hmax
        mesh_N = round((prob.xmax-prob.xmin)/prob.hmax);
        
        grid = setupUniformGrid(prob, mesh_N);
        % 
        % p = grid.p;
        % e = grid.e;
        % t = grid.t;
        % 
        % % ref=[];
        % ref=[10,10;5,5;2,2;1,1]';
        % % ref=[10,10;5,5]';
        % % ref=[5,5]';
        % for i=1:size(ref,2)
        %     tcenters = 1/3*(p(:,t(1,:))+p(:,t(2,:))+p(:,t(3,:)));
        %     % tidx = find(abs(tcenters(1,:)) <= ref(1,i) &abs(tcenters(2,:))<= ref(2,i) );
        %     tidx = find(sqrt(tcenters(1,:).^2+tcenters(2,:).^2) <= ref(1,i));
        %     t(4,:) = ones(1,size(t,2));
        %     [p,e,t]=refinemesh(prob.geom,p,e,t,tidx');
        % end
        % grid.p = p; grid.e = e; grid.t = t;
        % 
        % TODO This line was uncommented and seems to cause an error!
        % grid2 = setupIrregularGrid(prob, grid);
    
    
        % DEBUG
        % plotMesh(grid)
        % grid2 = setupIrregularGrid(prob, grid);
        % max(abs(grid2.muT-grid.muT))
        % max(abs(grid2.muOmega-grid.muOmega))
        % max(max(abs(grid2.muGamma-grid.muGamma)))
        % max(max(abs(grid2.D0-grid.D0)))
        % for i=1:2, max(max(abs(grid2.D1{i}-grid.D1{i}))), end
        % for i=1:2, max(max(abs(grid2.D1_up{i}-grid.D1_up{i}))), end
        % for i=1:4, max(max(abs(grid2.D2{i}-grid.D2{i}))), end
        % all(grid.inner==grid2.inner)
        % all(grid.bnd==grid2.bnd)
        % keyboard
    
    case 'irregular'
    
        grid = setupIrregularGrid(prob);
    	
    otherwise
    	'prob.mesh_type in setupGrid not known'
    end
    
    grid.dim = prob.dim;
    
    
    
    % -------------------------------------------
    % NOTES:
    % -------------------------------------------
    %
    % I'm wondering, wheather the matlab initmesh routine internally calls
    % the delaunayTriangulation routine.
    % It seems to be that 
    % dt.ConnectivityList ~ t(1:3,:)' (same triangles, maybe nodes permuted)
    % Until I know this, I will overwrite the triangulation generated by initmesh.
    % Use Delaunay triangulation
    
    
    % Introduce constraints
    % TODO: The constraints are not represented in the Voronoi diagram
    % c = [1,2;2,3;3,4;4,1];
    % dt = delaunayTriangulation(p(1,:)',p(2,:)',c);
    % Eventually dt.convexHull could help
else
    grid = setupIrregularGrid(prob,[], pnts, e, t);
    grid.dim = prob.dim;
end