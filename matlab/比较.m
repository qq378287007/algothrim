close all; clear; clc;

myfun = @(a,H) a(1) + (a(2) .* H + a(3)) ./ (H .* H + a(4) .* H + a(5)) ...
    + (a(6) .* H + a(7)) ./ (H .* H + a(8) .* H + a(9));

D = readmatrix('ChinaSteel_35CS250H.csv');
H = D(:, 1);
B = D(:, 2);

beta0 = ones(9, 1);

beta1 = nlinfit(H, B, myfun, beta0);
y1 = myfun(beta1, H);

beta2 = mynlinfit(H, B, myfun, beta0);
y2 = myfun(beta2, H);


figure();
plot(H, B, 'bo-');
plot(H, y1, 'r*-');
plot(H, y2, 'g*-');



title("BH Curve");
xlabel("H");
ylabel("B");
legend("Initial Data", "Fit Data", 'Location', 'SouthEast');


