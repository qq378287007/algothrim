close all; clear; clc;

myfun = @(beta, X) beta(1) + (beta(2) .* X ...
    + beta(3)) ./ (X .* X + beta(4) .* X + beta(5)) ...
    + (beta(6) .* X + beta(7)) ./ (X .* X + beta(8) .* X + beta(9));

file = 'D23_50.csv';
file = 'DW310_35.csv';
file = 'FLN8.csv';%bad
file = 'FLNG28.csv';%bad
file = 'M19_24G.csv';%bad

files = ls('*.csv');
for i = 1 : size(files, 1)
    file = files(i, :);
    file = strtrim(file);
    
    D = readmatrix(file);
    H = D(2:end, 2);
    B = D(2:end, 3);
    
    beta0 = ones(9, 1);
    
    beta1 = nlinfit(H, B, myfun, beta0);
    y1 = myfun(beta1, H);
    
    figure();
    hold on;
    plot(H, B, 'bo-');
    plot(H, y1, 'r*-');
    
    print(gcf, '-dpng', [file, '.png']);
    close;
end

