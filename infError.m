function [inf_error] = infError(v, vexact)

inf_error = max(abs(vexact - v));

end % function [inf_error] = infError(v, vexact)
