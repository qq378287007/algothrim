number = 20;
x = 10 * rand(1, number);
y = 10 * rand(1, number);

tri = delaunay(x, y);

hold on;
plot(x, y, 'r*');
for ii = 1 : size(tri, 1)
    plot([x(tri(ii, 1)) x(tri(ii, 2))], [y(tri(ii, 1)) y(tri(ii, 2))], 'b');
    plot([x(tri(ii, 2)) x(tri(ii, 3))], [y(tri(ii, 2)) y(tri(ii, 3))], 'b');
    plot([x(tri(ii, 3)) x(tri(ii, 1))], [y(tri(ii, 3)) y(tri(ii, 1))], 'b');
end
set(gca, 'box', 'on');
print(gcf, '-dpng', 'delaunary.png');



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
print(gcf, '-dpng', 'delaunary.png');








