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

