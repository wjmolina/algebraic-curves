function power_base_out = power_base(i, j, r, s, l)
[x_grid, y_grid] = meshgrid(linspace(0, l, s), linspace(0, l, r));
power_base_out = x_grid .^ i .* y_grid .^ j;
end