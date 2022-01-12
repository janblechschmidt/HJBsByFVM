function B = getFirstOrderMatrix(prob, grid, coef, upwinding, u)
helperFunctions;

    if nargin < 4
        upwinding = 0;
        u = [];
    elseif nargin < 5
        u = [];
    end

    if isstruct(coef)
        vals = evaluateCoefficient(prob, grid, coef.convection, u);
    else
        [m,n] = size(coef);
        for i=1:n
            vals{i} = coef(:,i);
        end
    end

    if upwinding
        B = sparse(grid.N, grid.N);
        for i = 1:grid.dim
            Bforw = spdiags(grid.muOmega.*posPart(vals{i}), [0], grid.N, grid.N) * grid.D1forw{i};
            Bback = spdiags(grid.muOmega.*negPart(vals{i}), [0], grid.N, grid.N) * grid.D1back{i};
            B = B + Bforw - Bback;
        end
        
        % This is the version that works but has not this nice property that the directions decouple
        % K_total = sparse(grid.N,grid.N);
        % for i = 1:grid.dim
        %     tmp = vals{i};
        %     K_total = K_total + spdiags(grid.muOmega.*tmp, [0], grid.N, grid.N) * grid.D1_up{i};
        % end
        % ind = (K_total < 0);
        % noind = (K_total > 0);

        % applyPartitioning = 0;
        % if ~applyPartitioning
        %     B = sparse(grid.N, grid.N);
        %     for i = 1:grid.dim
        %         tmp = spdiags(grid.muOmega.*vals{i}, [0], grid.N, grid.N) * grid.D1_up{i};
        %         B = B + tmp .* noind + spdiags(sum(tmp.*(ind),2), [0], grid.N, grid.N);
        %     end

        % else
        %     B = sparse(grid.N, grid.N);
        %     sections = ceil(grid.N/5000);

        %     partition = round(linspace(0,1,sections+1)*grid.N);
        %     for j = 1:sections
        %         fprintf('.')
        %         idx = (partition(j)+1):partition(j+1);
        %         for i = 1:grid.dim
        %             tmp = spdiags(grid.muOmega(idx).*vals{i}(idx), [0], length(idx), length(idx)) * grid.D1_up{i}(idx,:);
        %             B(idx,:) = B(idx,:) + tmp .* ~ind(idx,:) + spdiags(sum(tmp.*(ind(idx,:)),2), [partition(j)], length(idx), grid.N);
        %         end
        %     end
        %     fprintf('\n')

        % end



        % The following variant is even more memory consuming
        % B = sparse(grid.N, grid.N);
        % for i = 1:grid.dim
        %     T = spdiags(vals{i}, [0], grid.N, grid.N) * grid.D1_up{i};
        %     Toffdiag = T;
        %     Toffdiag(ind)=0;
        %     T(~ind)=0;
        %     T=spdiags(sum(T,2), [0], grid.N, grid.N);
        %     B = B + Toffdiag + T;
        % end
        % max(max(abs(B2-B)))
        % keyboard

    else
        B = sparse(grid.N, grid.N);
        for i = 1:grid.dim
            B = B + spdiags(grid.muOmega .* vals{i}, [0], grid.N, grid.N) * grid.D1{i};
        end

    end


end % function B = getFirstOrderMatrix(prob, grid, coef, upwinding, u)
