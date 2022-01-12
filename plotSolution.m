function plotSolution(prob, grid, coef, opt, sol)
if opt.plot_solution
    sfigure(3); clf;
	plotOnTriangulation(grid.p, grid.t' ,sol.v);
	% title('Value function');
    xlabel('$x^1$')
    ylabel('$x^2$')
    view([40 30])
    axis equal
    axis tight
    % fnamebase = 'MAT_CircleDomain';
    % fname = sprintf('./tex/%s_value_%d.png', fnamebase, grid.N);
    % print(gcf,fname,'-dpng','-r300'); 

    if isfield(sol,'u')
        for i = 1:prob.control_dim
            if isfield(coef, 'solution')
                sfigure(5); clf;
                subplot(prob.control_dim, 2, i);
                plotOnVoronoiGrid(grid, sol.u(:,i));
                title('Control function')
                subplot(prob.control_dim, 2, prob.control_dim+i);
                uexact = coef.optimal_control(grid.p');
                plotOnVoronoiGrid(grid, abs(sol.u(:,i)-uexact(:,i)));
                title('Control error')
            else
                sfigure(5+i); clf;
                % plot(prob.control_dim, 1, i);
                plotOnVoronoiGrid(grid, sol.u(:,i));
				% plotOnTriangulation(grid.p, grid.t' ,sol.u(:,i));
                % title('Control function')
                % view([32,32])
                view(2)
                colorbar('eastoutside')
                % view([140+180,40])
                xlabel('$x^1$')
                ylabel('$x^2$')

            end
            axis equal
            axis tight
            % caxis([-4, 4])
            % fname = sprintf('./tex/%s_control%d_%d.png', fnamebase, i, grid.N);
            % print(gcf,fname,'-dpng','-r300'); 
        end
    end % if isfield(sol,'u')

    if isfield(sol,'pi')
        sfigure(5); clf;
        plotOnVoronoiGrid(grid, sol.pi);
        title('Active set')
    end % if isfield(sol,'pi')

    if isfield(coef,'solution')
        sfigure(4); clf;
        title('Solution')
        subplot(2,1,1);
        if isfield(grid, 'time')
            vexact = coef.solution(grid.time,grid.p');
        else
            if isfield(prob,'Tmax')
                vexact = coef.solution(prob.Tmin, grid.p');
            else
                vexact = coef.solution(grid.p');
            end
            %TODO Diese ganzen Abfragen sollten sich ändern.

        end
        plotOnTriangulation(grid.p, grid.t', vexact);
        title('Exact solution');
        subplot(2,1,2);
        plotOnTriangulation(grid.p, grid.t', abs(vexact - sol.v));
        title('Error to exact solution');
        fprintf('Maximum error: %6.4e\n', max(abs(vexact - sol.v)));
    end
    fprintf('Plotted the solution. Hit a key to proceed...\n')
    pause
end
end % function plotSolution(prob, grid, coef, opt, sol)
