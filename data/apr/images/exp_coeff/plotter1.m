clear;
clc;
x_grid = - 15 : 15;
c_array = power_expansion_coefficients(6, - 12, 12);
% c_array = separable_expansion_coefficients(6, - 12, 12);
for k = 6 : - 1 : 0
    if k == 0
        m = 'o';
    elseif k == 1
        m = '*';
    elseif k == 2
        m = 'v';
    elseif k == 3
        m = 'x';
    elseif k == 4
        m = '^';
    elseif k == 5
        m = 's';
    elseif k == 6
        m = '>';
    end
    plot(x_grid, abs(c_array(1 + k, :)) + 1, 'Marker', m);
    set(gca, 'YScale', 'log');
    hold on;
end
legend('m = 6', 'm = 5', 'm = 4', 'm = 3', 'm = 2', 'm = 1', 'm = 0');
title('Power Expansion Coefficients \mu');
xlabel('k');
ylabel('abs(\mu_{k}^{m, 6}) + 1 (log scale)');