close all;
clear;
clc;

D = readmatrix('ChinaSteel_35CS250H.csv');
D = readmatrix('ChinaSteel_50CS350.csv');
D = readmatrix('DW415.csv');
D = readmatrix('JFE_Steel_35JN210.csv');
%D = readmatrix('JFE_Steel_35JN360.csv');
H = D(:, 1);
B = D(:, 2);
myfun = @(a,H) a(1) + (a(2) .* H + a(3)) ./ (H .* H + a(4) .* H + a(5)) ...
    + (a(6) .* H + a(7)) ./ (H .* H + a(8) .* H + a(9));

beta0 = ones(9, 1);
beta1 = nlinfit(H, B, myfun, beta0);
y = myfun(beta1, H);
rmse = sqrt(mean((B-y).^2));

figure();
plot(H, B, 'bo-', H, y, 'r*-');
title("BH Curve");
xlabel("H");
ylabel("B");
legend("Initial Data", "Fit Data", 'Location', 'SouthEast');


%
beta2 =  [2.19028 ...
     22.3963  ...
    -2911.36  ...
    -96.5255  ...
     5020.67  ...
    -8289.45  ...
-1.7912e+006  ...
     12291.3  ...
1.10915e+006];

y = myfun(beta2, H);
rmse = sqrt(mean((B-y).^2));

figure();
plot(H, B, 'bo-', H, y, 'r*-');
title("BH Curve");
xlabel("H");
ylabel("B");
legend("Initial Data", "Fit Data", 'Location', 'SouthEast');
