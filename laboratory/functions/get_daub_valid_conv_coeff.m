function [c_array, f_array] = get_daub_valid_conv_coeff(n, a, b)
    n = n + 1;
    if mod(n, 2) == 0
        a = a - 1;
    end
    b = b + 1;
    p_array = 2 * dbwavf(['db', num2str(n)]);
    for k = 0 : n - 1
        c_array(k + 1, :) = pol_coeffs(n, k, p_array, a, b);
    end
    f_array = scaling_integers(p_array);
end