function [V, C, inside] = computeBoundedVoronoiCells(prob, grid, V, C, DV)

method = 1;

if method==1 % New method - works for polygonal bounded domains
    
    inside = zeros(grid.N,1);
    for i = 1:grid.N
    	cell = C{i};
    	if ~ismember(1,cell) % The index 1 is V(1,:) = [Inf, Inf]
    		inside(i) = 1;
            C{i} = C{i}-1;
        end
    end
    % Only vertices that belong to boundary edges have unbounded Voronoi cells
    % x_idx = unique(e(:));
    
    nv = size(V,1);

    % Get number of boundary edges
    ne = size(prob.geom,2);

    % Initialize counter for boundary values
    nev = [];

    % Add new boundary points to vertex list
    for i = 1:ne
        ei = grid.e(5,:)==i;
        e1 = grid.e(1,ei);
        e2 = grid.e(2,ei);

        % Compute new vertices along boundary
        Vnew = 0.5*(grid.p(:,e1)+grid.p(:,e2));

        % Append them to vertex list
        V = [V;Vnew'];

        e_inner = e1(2:end);
        ne_inner = length(e_inner);

        for j = 1:ne_inner
            k = e_inner(j);
            idx = nv+sum(nev)+j;
            [~,ord] = min(sum((V(idx,:)'-V(C{k}(2:end),:)').^2,1));
            if ord == 1
                idx = [idx,C{k}(2:+1:end)];
            else
                idx = [idx,C{k}(end:-1:2)];
            end
            idx = [idx, nv+sum(nev)+j+1];

            C{k} = idx-1; 
        end

        % Store number of new elements
        nev = [nev,size(Vnew,2)];

    end
    nev_total = sum(nev);

    % Add corners to vertex list
    p1 = find(grid.e(3,:)==0);
    p2 = find(grid.e(4,:)==1);
    Vnew = grid.p(:,grid.e(1,p1));
    V = [V;Vnew'];
    
    % Loop through corners
    for k = 1:ne
        pi = p1(k);
        i = grid.e(1,pi);
        e1 = grid.e(5,pi);
        

        % Start with corner
        idx = [nv+nev_total+k];

        % Next comes one new vertex along the edge belonging to pi
        new_idx = nv+sum(nev(1:e1-1))+1;
        idx = [idx,new_idx];

        % From this point, we look for the closest point in V(C{i},:)
        [~,ord] = min(sum((V(new_idx,:)'-V(C{i}(2:end),:)').^2,1));
        if ord == 1
            idx = [idx,C{i}(2:+1:end)];
        else
            idx = [idx,C{i}(end:-1:2)];
        end

        % Find intersection with other boundary
        e2 = grid.e(5,p2(grid.e(2,p2)==i));
        new_idx = nv+sum(nev(1:e2));
        idx = [idx, new_idx];

        % Set vertex indices to cell array
        % Offset of -1 is necessary because first variable is (inf,inf) and is removed
        C{i} = idx-1; 

        % Check result
        V(C{i}+1,:);
    end

    V = V(2:end,:);
    % keyboard

    % e_idx = unique(grid.e(1:2,:));
    % c_idx = grid.e(1,p1);
    % for i=e_idx
    %     if isempty(intersect(i,c_idx)) % Not a corner node
    %         i
    %     end
    % end
    % keyboard
    % 
    % Dx = computeDistanceToBoundary()
    % 
    % inside = zeros(grid.N,1);
    % m = size(V,1);
    % p = grid.p;
    % for i = 1:grid.N
    % 	cell = C{i}
    % 	
    %     if ~ismember(1,cell) % The index 1 is V(1,:) = [Inf, Inf]
    % 		inside(i) = 1;
    % 
    %     else
    %         xi = p(:,i)
    %         V(cell,:)
    %         keyboard
    % 
    %     end
    % 
    % 
    % end
    % 

else % Previous method - works only for rectangular domains
    ABS_ERROR = 1e-10;
    
    inside = zeros(grid.N,1);
    m = size(V,1);
    p = grid.p;
    for i = 1:grid.N
    	cell = C{i};
    
    	if ~ismember(1,cell) % The index 1 is V(1,:) = [Inf, Inf]
    		inside(i) = 1;
    
    	else
    
    		xi = p(:,i)';
    		% Corners have to be handled with care
    		corners = [prob.xmin, prob.xmin, prob.xmax, prob.xmax;
    							 prob.ymin, prob.ymax, prob.ymin, prob.ymax];
    		corner_index = find(all(xi'==corners));
    
    		if ~isempty(corner_index)
    			wx = V(cell(2),:);
    			valx = wx;
    			if any(corner_index == [1,4])
    				valx(2) = corners(2,corner_index);
    			else
    				valx(1) = corners(1,corner_index);
    			end
    
    			wy = V(cell(end),:);
    			valy = wy;
    			
    			if any(corner_index == [1,4])
    				valy(1) = corners(1,corner_index);
    			else
    				valy(2) = corners(2,corner_index);
    			end
    
    			valc = corners(:,corner_index)';
    
    			% Check, whether points already belong to V
    			wx_idx = find(sum((V - valx).^2,2) < ABS_ERROR);
    			wy_idx = find(sum((V - valy).^2,2) < ABS_ERROR);
    			wc_idx = find(sum((V - valc).^2,2) < ABS_ERROR);
    			
    			if isempty(wx_idx)
    				wx_idx = m+1;
    				V(wx_idx,:) = valx;
    				m = m+1;
    			end
    			C{i}(1) = wx_idx;
    
    			if isempty(wy_idx)
    				wy_idx = m+1;
    				V(wy_idx,:) = valy;
    				m = m+1;
    			end
    			C{i}(end+1) = wy_idx;
    
    			if isempty(wc_idx)
    				wc_idx = m+1;
    				V(wc_idx,:) = valc;
    				m = m+1;
    			end
    			C{i}(end+1) = wc_idx;
    
    		else
    			% Determine boundary, at which point xi lies
    			[~,bnd_idx] = determineClosestBndValue(prob, xi);
    			wa = V(cell(2),:);
    			we = V(cell(end),:);
    			switch bnd_idx
    			case 1
    				wa(1) = prob.xmin;
    				we(1) = prob.xmin;
    			case 2
    				wa(1) = prob.xmax;
    				we(1) = prob.xmax;
    			case 3
    				wa(2) = prob.ymin;
    				we(2) = prob.ymin;
    			case 4
    				wa(2) = prob.ymax;
    				we(2) = prob.ymax;
    			end
    
    			% Check, whether points already belong to V
    			wa_idx = find(sum((V - wa).^2,2) < ABS_ERROR);
    			we_idx = find(sum((V - we).^2,2) < ABS_ERROR);
    
    			if isempty(wa_idx)
    				wa_idx = m+1;
    				V(wa_idx,:) = wa;
    				m = m+1;
    			end
    			C{i}(1) = wa_idx;
    
    			if isempty(we_idx)
    				we_idx = m+1;
    				V(we_idx,:) = we;
    				m = m+1;
    			end
    			C{i}(end+1) = we_idx;
    		end
    	end
    end
end





end % function computeBoundedVoronoiCells(grid, V, C)
