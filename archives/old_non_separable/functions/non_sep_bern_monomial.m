function nsbm = non_sep_bern_monomial(domain_length, a, b, n, i, j)
[x_grid, y_grid] = meshgrid(linspace(0, domain_length, b), linspace(0, domain_length, a));
nsbm = zeros(a, b);
for p = 0 : n - i - j
    for q = 0 : n - i - j - p
        nsbm = nsbm + nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(p) * factorial(q) * factorial(n - i - j - p - q)) * domain_length ^ (p - n) * (- 1) ^ (n - i - j - p) * x_grid .^ (q + i) .* y_grid .^ (n - i - j - p - q + j);
    end
end
end