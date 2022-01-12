function [geom,gd] = getDomain(domainLabel, xmin, xmax, ymin, ymax)
if strcmpi(domainLabel, 'UnitSquare')
    gd = [3, 4, 0, 1, 1, 0, 0, 0, 1, 1]';
    gm = [gd];
    sf = 'SQ1';
    ns = char('SQ1');
    ns = ns';
    geom = decsg(gm, sf, ns);
end
if strcmpi(domainLabel, 'Rectangle')
    gd = [3, 4, xmin, xmax, xmax, xmin, ymin, ymin, ymax, ymax]';
    gm = [gd];
    sf = 'R1';
    ns = char('R1');
    ns = ns';
    geom = decsg(gm, sf, ns);
end
if strcmpi(domainLabel, 'Diamond')
    xmid = 0.5*(xmax+xmin);
    ymid = 0.5*(ymax+ymin);
    gd = [2, 4, xmin, xmid, xmax, xmid, ymid, ymin, ymid, ymax]';
    gm = [gd];
    sf = 'P1';
    ns = char('P1');
    ns = ns';
    geom = decsg(gm, sf, ns);
end

end
