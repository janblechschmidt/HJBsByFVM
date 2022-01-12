function c2f = setupFineToCoarseMatrix(gfine, g)

n = g.N;
nfine = gfine.N;
switch g.dim
case 1
    row = 1:n;
    col = 1:2:nfine;
    f2c = sparse(row(:), col(:), ones(n,1), n, nfine);

    row = [1:n-1,2:n];
    col = [2:2:nfine-1, 2:2:nfine-1];
    f2c = f2c + sparse(row(:), col(:), 0.5*ones(2*(n-1),1), n, nfine);
case 2
    
    % Define numbers of nodes as they are in the mesh
    I = reshape(1:n, sqrt(n), sqrt(n));
    Ifine = reshape(1:nfine, sqrt(nfine), sqrt(nfine));
    
    % Identify the previous nodes with the new ones
    col = Ifine(1:2:end,1:2:end);
    row = I;
    f2c = sparse(row(:), col(:), ones(n,1), n, nfine);
    
    % Find the nodes that are interpolated by values above and below
    col = Ifine(2:2:end,1:2:end); 
    row1 = I(1:end-1,:);
    row2 = I(2:end,:);
    f2c = f2c + sparse([row1(:);row2(:)],[col(:);col(:)],0.5*ones(2*(n-sqrt(n)),1),n, nfine);
    
    % Find the nodes that are interpolated by values left and right
    col = Ifine(1:2:end,2:2:end);
    row1 = I(:,1:end-1);
    row2 = I(:,2:end);
    
    f2c = f2c + sparse([row1(:);row2(:)],[col(:);col(:)],0.5*ones(2*(n-sqrt(n)),1),n, nfine);
    
    % Find the nodes that are interpolated by diagonal values
    col = Ifine(2:2:end,2:2:end);
    row1 = I(1:end-1,1:end-1);
    row2 = I(2:end,2:end);
    
    f2c = f2c + sparse([row1(:);row2(:)],[col(:);col(:)],0.5*ones(2*((sqrt(n)-1)^2),1),n, nfine);
    
end

c2f = f2c';

end % function f2c = setupFineToCoarseMatrix(g, gfine)
