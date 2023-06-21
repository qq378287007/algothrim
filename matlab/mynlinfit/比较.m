close all; clear; clc;

myfun = @(beta, X) beta(1) + (beta(2) .* X ... 
    + beta(3)) ./ (X .* X + beta(4) .* X + beta(5)) ...
    + (beta(6) .* X + beta(7)) ./ (X .* X + beta(8) .* X + beta(9));

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

figure();
plot(H, B, 'bo-');
plot(H, y2, 'g*-');

%{
title("BH Curve");
xlabel("H");
ylabel("B");
legend("Initial Data", "Fit Data", 'Location', 'SouthEast');
%}

all(y1 == y2)

X = H;
Y = B;
model = myfun;
beta = beta0;
maxiter = 200;
