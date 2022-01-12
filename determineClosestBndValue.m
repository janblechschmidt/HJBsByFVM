function [val,arg] = closestBndValue(prob, w)

toxmin = norm(w(1) - prob.xmin);
toxmax = norm(w(1) - prob.xmax);
toymin = norm(w(2) - prob.ymin);
toymax = norm(w(2) - prob.ymax);
[minval,arg] = min([toxmin, toxmax, toymin, toymax]);
arg = find(minval==[toxmin, toxmax, toymin, toymax]);
j = 1;
for i=arg
	val{j} = w;
	switch i
	case 1
		val{j}(1) = prob.xmin;
	case 2
		val{j}(1) = prob.xmax;
	case 3
		val{j}(2) = prob.ymin;
	case 4
		val{j}(2) = prob.ymax;
	end
	j = j + 1;
end

end
