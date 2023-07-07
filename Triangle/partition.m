fileID = fopen('out2.txt');
p_num = fscanf(fileID, '%d', 1);
points = fscanf(fileID, '%f', [2 p_num]);
points = points';

edge_num = fscanf(fileID, '%d', 1);
edge = fscanf(fileID, '%f', [4 edge_num]);
edge = edge';
fclose(fileID);


hold on;
plot(points(:, 1), points(:, 2), 'r*');
for ii = 1 : size(edge, 1)
    plot([edge(ii, 1) edge(ii, 3)], [edge(ii, 2) edge(ii, 4)], 'b');
end
set(gca, 'box', 'on');
print(gcf, '-dpng', 'partition.png');

