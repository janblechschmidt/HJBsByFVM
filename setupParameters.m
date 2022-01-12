% This file sets all parameters for the energy storage problem
function p = setupParameters(p)

if p.domain == 'UnitSquare'
        p.geom = [  2     2     2     2;
        0     1     1     0;
        1     1     0     0;
        0     0     1     1;
        0     1     1     0;
        1     1     1     1;
        0     0     0     0];

end

if strcmp( p.mode, 'MinimumArrivalTime' )
    p = parseFEniCSxml('HJBMinimumArrivalTime.xml',p);

    if p.dim == 2
        p.geom = [  2     2     2     2;
        0     1     1     0;
        1     1     0     0;
        0     0     1     1;
        0     1     1     0;
        1     1     1     1;
        0     0     0     0];
    end %if p.dim == 2

end % if strcmp( p.mode, 'MinimumArrivalTime' )

if strcmp( p.mode, 'MinimumArrivalTime_xcontrol' )
    p = parseFEniCSxml('HJBMinimumArrivalTime.xml',p);

    if p.dim == 2
        p.geom = [  2     2     2     2;
        0     1     1     0;
        1     1     0     0;
        0     0     1     1;
        0     1     1     0;
        1     1     1     1;
        0     0     0     0];
    end %if p.dim == 2

end % if strcmp( p.mode, 'MinimumArrivalTime_xcontrol' )


end % function p = setupParameters(p)

% Box with hole domain
%  p.geom = [2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000;
%  1.5000    1.5000    0.5000    0.5000   -0.5000   -1.5000   -1.5000   -0.5000;
%  1.5000   -1.5000    0.5000   -0.5000    0.5000   -1.5000    1.5000   -0.5000;
%  1.0000   -1.0000    0.4000   -0.4000    0.4000   -1.0000    1.0000   -0.4000;
% -1.0000   -1.0000   -0.4000   -0.4000    0.4000    1.0000    1.0000    0.4000;
%       0         0    1.0000    1.0000    1.0000         0         0    1.0000;
%  1.0000    1.0000         0         0         0    1.0000    1.0000         0];
