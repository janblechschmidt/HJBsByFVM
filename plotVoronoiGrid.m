function ax = plotVoronoiGrid(grid, ax)

% title('Corresponding bounded dual grid')
xlabel('$x$'), ylabel('$y$')
for i = 1:grid.N
	cell = grid.C{i};
	plot(grid.V([cell,cell(1)],1),grid.V([cell,cell(1)],2),'k-')
end

if grid.N < 100
	plot(grid.p(1,:)',grid.p(2,:)','b.')
	text(grid.p(1,:)',grid.p(2,:)',num2str((1:grid.N)'))
end

if nargin < 2
	widen_axes(0.1);
	ax = axis();
else
	axis(ax);
end
axis equal
axis tight

end % function ax = plotVoronoiGrid(grid, ax)
