function power_coeff_table = power_expansion_coefficients(n, a, b)
    n_total = n + 1;
    p_array = zeros(1, n_total + 1);
    for i = 0 : n_total
        p_array(1 + i) = nchoosek(n_total, i);
    end
    p_array = p_array * 2 ^ (- n_total + 1);
    power_coeff_table = pol_coeffs(n_total, 0, p_array, a, b);
    for k = 1 : n
        power_coeff_table(1 + k, :) = pol_coeffs(n_total, k, p_array, a, b);
    end
end