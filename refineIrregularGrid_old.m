function grid = refineIrregularGrid(prob, grid)
% This function has to update all grid-relevant quantities, both for the primal and dual mesh
e = grid.edges;
p = grid.p;
t = grid.t(1:3,:);
np = size(grid.p,2)
% x coordinates of new points
px = 0.5*(p(1,e(1,:))+p(1,e(2,:)));
py = 0.5*(p(2,e(1,:))+p(2,e(2,:)));
new_np = size(e,2);
new.p = [grid.p,[px;py]];

% Compute new edges (they are not directed, but always have the lower index in the first row)
new.edges = reshape([grid.edges(1,:); (1:new_np)+np; grid.edges(2,:); (1:new_np)+np],2,[]);
t_new = [];
keyboard
for i = 1:grid.N_tri
    ti = t(:,i);

    e_ti_1 = sort(ti([1,2]));
    e_ti_2 = sort(ti([2,3]));
    e_ti_3 = sort(ti([3,1]));

    idx_1 = find(all(grid.edges==e_ti_1))+np;
    idx_2 = find(all(grid.edges==e_ti_2))+np;
    idx_3 = find(all(grid.edges==e_ti_3))+np;
    
    % Append new triangles with the same orientation as the old ones
    t_new = [t_new, [ti(1);idx_1;idx_3], [ti(2);idx_2;idx_1], [ti(3);idx_3;idx_2], [idx_1;idx_2;idx_3] ];
    
end
new.t = t_new

volT = grid.muT*[.25,.25,.25,.25];
new.muT = reshape(volT',[],1);

% keyboard

end % function refineIrregularGrid(prob, grid)
