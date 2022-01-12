% Updated 2021-11-08
% Load helper functions
clear all, clc
helperFunctions

prob.dim = 2;
prob.mode = 'MinimumArrivalTime';
% prob.mode = 'CinftySolTime';
prob.outputdir = 'output2d';
prob.trivialControl = 0;
prob.domain = 'UnitSquare';
% prob.mesh_type = 'irregular';
prob.mesh_type = 'regular';
% prob.dual_mesh = 'voronoi';
% prob.grid_type = 'left'
% prob.grid_type = 'right';

prob.DirichletBoundary = @(onBoundary,x) onBoundary;
prob.hmax = 0.5; % Coarse initial grid
% prob.hmax = 0.25; % Coarse initial grid
% prob.hmax = 0.1;
% prob.hmax = 0.13; % Fine initial grid
% prob.hmax = 0.09; % Fine initial grid
% prob.hmax = 0.03; % Fine initial grid

prob.xmin = 0;
prob.xmax = 1;
prob.ymin = 0;
prob.ymax = 1;

prob = setupParameters(prob);
coef = setupCoefficients(prob);
grid = setupGrid(prob);
% Required input:
% p - Points in mesh
% e - Edges in mesh
% t - Triangles in mesh


% Here begins the export of the primal grid
fid = fopen('./tex/primalGrid.tex','w');

fprintf(fid,'\\begin{tikzpicture}[scale=0.65, myLine/.style={blue}]');
fprintf(fid,'\n\\begin{axis}[xmin=-.1,   xmax=1.1, ymin=-.1,   ymax=1.1]');
for i = 1:size(grid.t,2)

    fprintf(fid,'\n\\addplot[myLine] coordinates {\n');
    for j = [1,2,3,1]
        fprintf(fid,'(%8.6f, %8.6f)\n', grid.p(:,grid.t(j,i)) );
    end
    fprintf(fid,'};');
end
fprintf(fid,'\n\\end{axis}');
fprintf(fid,'\n\\end{tikzpicture}\n');
fclose(fid);


% Here begins the export of the dual grid
fid = fopen('./tex/dualVoronoiGrid.tex','w');

fprintf(fid,'\\begin{tikzpicture}[scale=0.65, myLine/.style={blue}, primalMesh/.style={gray, opacity=0.3}]');
fprintf(fid,'\n\\begin{axis}[xmin=-.1,   xmax=1.1, ymin=-.1,   ymax=1.1]');
for i = 1:size(grid.C,2)

    fprintf(fid,'\n\\addplot[myLine] coordinates {\n');
    for j = [grid.C{i}, grid.C{i}(1)]
        fprintf(fid,'(%8.6f, %8.6f)\n', grid.V(j,:) );
    end
    fprintf(fid,'};');
end
for i = 1:size(grid.t,2)

    fprintf(fid,'\n\\addplot[primalMesh] coordinates {\n');
    for j = [1,2,3,1]
        fprintf(fid,'(%8.6f, %8.6f)\n', grid.p(:,grid.t(j,i)) );
    end
    fprintf(fid,'};');
end
fprintf(fid,'\n\\end{axis}');
fprintf(fid,'\n\\end{tikzpicture}\n');
fclose(fid);

prob.dual_mesh = 'centroid';
prob = setupParameters(prob);
coef = setupCoefficients(prob);
grid = setupGrid(prob);

% Here begins the export of the dual grid
fid = fopen('./tex/dualCentroidGrid.tex','w');

fprintf(fid,'\\begin{tikzpicture}[scale=0.65, myLine/.style={blue}, primalMesh/.style={gray, opacity=0.3}]');
fprintf(fid,'\n\\begin{axis}[xmin=-.1,   xmax=1.1, ymin=-.1,   ymax=1.1]');
for i = 1:size(grid.C,2)

    fprintf(fid,'\n\\addplot[myLine] coordinates {\n');
    for j = [grid.C{i}, grid.C{i}(1)]
        fprintf(fid,'(%8.6f, %8.6f)\n', grid.V(j,:) );
    end
    fprintf(fid,'};');
end
for i = 1:size(grid.t,2)

    fprintf(fid,'\n\\addplot[primalMesh] coordinates {\n');
    for j = [1,2,3,1]
        fprintf(fid,'(%8.6f, %8.6f)\n', grid.p(:,grid.t(j,i)) );
    end
    fprintf(fid,'};');
end
fprintf(fid,'\n\\end{axis}');
fprintf(fid,'\n\\end{tikzpicture}\n');
fclose(fid);
