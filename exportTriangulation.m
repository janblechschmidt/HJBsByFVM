function exportTriangulation(fname, grid)
% Here begins the export of the primal grid
fid = fopen(fname,'w');

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

end
