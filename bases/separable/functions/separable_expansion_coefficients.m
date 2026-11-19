function sec = separable_expansion_coefficients(n, a, b)
    pec = power_expansion_coefficients(n, a, b);
    sec = zeros(1 + n, size(pec, 2));
    for k = 0 : n
        for r = 0 : k
            for s = 0 : n - k
                sec(1 + k, :) = sec(1 + k, :) + nchoosek(n, k) * nchoosek(k, r) * nchoosek(n - k, s) * (b - a) ^ (- n) * a ^ (k - r) * b ^ (n - k - s) * (- 1) ^ (k + r + s) * pec(1 + r + s, :);
            end
        end
    end
end