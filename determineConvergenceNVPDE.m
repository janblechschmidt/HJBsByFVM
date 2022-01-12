function [sol, conv, grid] = determineConvergenceNVPDE(prob, grid, coef, opt)

markedTriangles = [];
markedPoints = [];

% Initialize lists for errors
conv.N = [];
conv.hmax = [];
conv.inf_error = [];
conv.L2_error = [];

% Initialize refinement counter
num_ref = 1;
ref_done = 0;

while ~ref_done
	fprintf('\n\nComputation on level %d\n\n', num_ref);
	
    % TODO: Implement own routine for refinement
    switch opt.refinementType
    case 'global' % Global refinement
        if num_ref > 1
            grid = refineGrid(prob, grid, [], []);
            if strfind(prob.class, 'parabolic')
                prob.dtmax = prob.dtmax/2;
                fprintf('Reduced dtmax to %6.4f\n', prob.dtmax);
            end
        end
    case 'local'  % Adaptive refinement
        if num_ref > 1
            
            % Estimate errors on triangles
            v = V{num_ref-1};
            % markTriangles(v, grid)
            keyboard
            error('Adaptive refinement has to be implemented\n')
        end
    end % switch opt.refinementType
   
    plotMesh(grid, opt);
   
    if strfind(prob.class, 'elliptic')
        if strfind(prob.class, 'HJB')
            keyboard
            sol = solveEllipticHJB(prob, grid, coef, opt);
        else
            sol = solveNVPDE(prob, grid, coef, opt);
        end
    else
        sol = solveParabolicProblem(prob, grid, coef, opt);
    end

    % Store grids (this is only necessary if there is no explicit solution)
    Grid{num_ref} = grid;

    % Deconstruct sol structure
    V{num_ref} = sol.v;
    if strfind(prob.class, 'HJB')
        U{num_ref} = sol.u;
    end

    % We plot the results, if opt.plot_solution == True 
    plotSolution(prob, grid, coef, opt, sol);

    if isfield(coef,'solution')
        if isfield(prob, 'Tmax')
            vexact = coef.solution(prob.Tmin, grid.p');
        else
            vexact = coef.solution(grid.p');
        end

        hmax = max(max(grid.muGamma));
        conv.N = [conv.N, grid.N];
        conv.hmax = [conv.hmax, hmax];
        conv.inf_error = [conv.inf_error, infError(sol.v, vexact)];
        conv.L2_error = [conv.L2_error, L2Error(grid, sol.v, vexact)];
    end
        
	switch opt.refinementType
	case 'no' % no refinement
		
		ref_done = 1;

	case 'global' % uniform refinement

		if num_ref > opt.number_of_refinements
			ref_done = 1;
        else
            num_ref = num_ref + 1;
		end

    case 'local'

        fprintf('Adaptive refinement has to be implemented\n');

		num_ref = num_ref + 1;
		if num_ref > opt.number_of_refinements
			ref_done = 1;
		end
    end % switch opt.refinement_type
end % while ~ref_done

% Handle the case without explicit solution
if ~isfield(coef,'solution')
    % Take the last solution as reference
    v_fine = sol.v;
    grid_fine = grid;
    for i = 1:num_ref-1
        grid = Grid{i};
        v = V{i};
        for j=i+1:num_ref
            v = Grid{j}.c2f*v;
        end
        hmax = max(max(grid.muGamma));
        conv.N = [conv.N, grid.N];
        conv.hmax = [conv.hmax, hmax];
        conv.inf_error = [conv.inf_error, infError(v, v_fine)];
        conv.L2_error = [conv.L2_error, L2Error(grid_fine, v, v_fine)];
    end
end

conv.inf_rate = ( log(conv.inf_error(2:end)./conv.inf_error(1:end-1)) ) ./ ( log(conv.hmax(2:end)./conv.hmax(1:end-1)) );
conv.L2_rate = ( log(conv.L2_error(2:end)./conv.L2_error(1:end-1)) ) ./ ( log(conv.hmax(2:end)./conv.hmax(1:end-1)) );


end % function conv = determineConvergenceNVPDE(prob, grid, coef, opt)
