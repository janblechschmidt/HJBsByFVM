function ax = plotTriangulation(grid, ax)
% p = dt.Points';
if nargin < 2
    figure(1); clf;
    ax = gca
end
p = grid.p;
t = grid.t(1:3,:);
N = size(p,2);
N_tri = size(t,2);
triplot(t', p(1,:), p(2,:))
% title('Triangulation')
xlabel('$x$'), ylabel('$y$')
hold on
if N < 100
    plot(p(1,:)',p(2,:)','b.')
    text(p(1,:)',p(2,:)',num2str((1:N)'))
    xc = 1/3*(p(1,t(1,:))+p(1,t(2,:))+p(1,t(3,:)));
    yc = 1/3*(p(2,t(1,:))+p(2,t(2,:))+p(2,t(3,:)));
    text(xc, yc, num2str((1:N_tri)'), 'Color','red')
end

if nargin < 2
	widen_axes(0.1);
	ax = axis();
else
	axis(ax);
end
axis equal
axis tight

hold off
