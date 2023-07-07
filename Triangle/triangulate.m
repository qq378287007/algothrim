fileID = fopen('out.txt');
p_num = fscanf(fileID, '%d', 1);
points = fscanf(fileID, '%f', [2 p_num]);
points = points';

tri_num = fscanf(fileID, '%d', 1);
tri = fscanf(fileID, '%f', [6 tri_num]);
tri = tri';
fclose(fileID);


hold on;
plot(points(:, 1), points(:, 2), 'r*');
for ii = 1 : size(tri, 1)
    plot([tri(ii, 1) tri(ii, 3)], [tri(ii, 2) tri(ii, 4)], 'b');
    plot([tri(ii, 3) tri(ii, 5)], [tri(ii, 4) tri(ii, 6)], 'b');
    plot([tri(ii, 5) tri(ii, 1)], [tri(ii, 6) tri(ii, 2)], 'b');
end
set(gca, 'box', 'on');
print(gcf, '-dpng', 'triangle.png');


