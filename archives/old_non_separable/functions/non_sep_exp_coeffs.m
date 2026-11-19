function nsec = non_sep_exp_coeffs(n, r, s)
c1 = dual_expansion_coefficients(n, r, 1);
c2 = dual_expansion_coefficients(n, s, 1);
nsec = zeros(n + 1, n + 1, size(c1, 2), size(c2, 2));
for i = 0 : n
    for j = 0 : n - i
        for a = 0 : n - i - j
            for b = 0 : n - i - j - a
                nsec(1 + i, 1 + j, :, :) = nsec(1 + i, 1 + j, :, :) + reshape(c1(1 + b + j, :)' * nchoosek(n, i) * nchoosek(n - i, j) * factorial(n - i - j) / (factorial(n - i - j - a - b) * factorial(a) * factorial(b)) * (- 1) ^ (a + b) * c2(1 + a + i, :), 1, 1, size(c1, 2), size(c2, 2));
            end
        end
    end
end
end