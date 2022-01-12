function [L2_error] = L2Error(grid, v, vexact)

L2_error = sqrt(sum(grid.muOmega.*((v-vexact).^2)));

end % function [inf_error] = infNorm(v, vexact)
