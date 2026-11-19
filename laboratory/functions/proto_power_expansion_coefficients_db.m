function c_array = proto_power_expansion_coefficients_db(n, a, b)
    p_array = 2 * dbwavf("db" + num2str(n));
    for k = 0 : n - 1
        c_array(1 + k, :) = pol_coeffs(n, k, p_array, a - 1, b + 1 - mod(n, 2));
    end
end