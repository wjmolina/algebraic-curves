function nsb = non_sep_bern(domain_length, a, b, n, i, j)
[x_grid, y_grid] = meshgrid(linspace(0, domain_length, b), linspace(0, domain_length, a));
nsb = domain_length ^ (- n) * nchoosek(n, i) * nchoosek(n - i , j) * x_grid .^ i .* y_grid .^ j .* (domain_length - x_grid - y_grid) .^ (n - i - j);
end