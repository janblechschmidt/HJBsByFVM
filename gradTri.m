function grad = gradTri(coords,v_local)
% Compute gradient of linear function on triangle
% with   x = coords(:,1),   y = coords(:,2),   z = coords(:,3),
% and  v(x) = v_local(1), v(y) = v_local(2), v(z) = v_local(3),

	dx1 = coords(1,2) - coords(1,1);
	dx2 = coords(1,3) - coords(1,1);
	dy1 = coords(2,2) - coords(2,1);
	dy2 = coords(2,3) - coords(2,1);
	dv1 = v_local(2) - v_local(1);
	dv2 = v_local(3) - v_local(1);

	grad = zeros(2,1);
	if dx1 ~= 0
		grad(2) = (dv2 - dv1 * dx2 / dx1) / (dy2 - dy1 * dx2 / dx1);
		grad(1) = dv1 / dx1 - grad(2) * dy1 / dx1;
	else
		if dx2 ~= 0
			grad(2) = (dv1 - dv2 * dx1 / dx2) / (dy1 - dy2 * dx1 / dx2);
			grad(1) = dv2 / dx2 - grad(2) * dy2 / dx2;
		else
			error('implement grad')
		end
	end

end % function grad = gradTri(coords,v_local)
