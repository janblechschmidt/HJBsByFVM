function grid = setupGridBoundaries(grid, prob)

    if isfield(prob, 'xmin')
        grid.lowbnd{1} = grid.p(1,:) <= prob.xmin + 1e-10;
        grid.upbnd{1} = grid.p(1,:) >= prob.xmax - 1e-10;
        grid.lowbnd{2} = grid.p(2,:) <= prob.ymin + 1e-10;
        grid.upbnd{2} = grid.p(2,:) >= prob.ymax - 1e-10;
    end
    
    x_s = grid.p(:,grid.e(1,:));
    x_e = grid.p(:,grid.e(2,:));
    Ndomain = flipud(x_s - x_e).*(-1).^grid.e(6:7,:);
    Ndomain = Ndomain./sqrt(sum(Ndomain.^2,1));
    Ndomain(3,:) = grid.e(1,:);
    Ndomain(4,:) = sqrt(sum((x_s-x_e).^2,1))/2;
    ndom = size(Ndomain,2);
    Ndomain = repmat(Ndomain,1,2);
    Ndomain(3,ndom+1:end) = grid.e(2,:);
    grid.n_boundary = Ndomain;

end
