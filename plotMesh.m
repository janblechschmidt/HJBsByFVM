function plotMesh(grid, opt)
if nargin < 2 | opt.plot_mesh 
    sfigure(1); clf, hold on
    ax = plotTriangulation(grid);
    % print(gcf,'./tex/WorstAssetTriangulation.png','-dpng','-r300'); 

	sfigure(2); clf, hold on
	plotVoronoiGrid(grid, ax);
    % print(gcf,'./tex/WorstAssetDualMesh.png','-dpng','-r300'); 
    fprintf('Plotted the meshes. Hit a key to proceed...\n')
    pause
end
end % function plotMesh(grid, opt)
