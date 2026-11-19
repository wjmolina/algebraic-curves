function deco = dual_expansion_coefficients(n, r, l)
    n = n + 1;
    r = r - 1;
    p_array = zeros(1, n + 1);
    for i = 0 : n
        p_array(1 + i) = nchoosek(n, i);
    end
    p_array = p_array * 2 ^ (- n + 1);
    deco = zeros(n, r + n - (mod(n, 2) == 0));
    for k = 0 : n - 1
        deco(1 + k, :) = pol_coeffs(n, k, p_array, 0, r) * (r / l) ^ (- k);
    end
end