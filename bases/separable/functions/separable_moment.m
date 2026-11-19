function moment_out = separable_moment(image_data, n, m, k, l, a, b, c, d)
    if nargin == 5
        a = 0;
        b = size(image_data, 2) - 1;
        c = 0;
        d = size(image_data, 1) - 1;
    end
    if k >= 0 && l >= 0 && n >= k && m >= l
        k1 = separable_expansion_coefficients(n, c, d);
        k2 = separable_expansion_coefficients(m, a, b);
        k1 = k1(k + 1, 1 + floor(n / 2) : end - floor(n / 2));
        k2 = k2(l + 1, 1 + floor(m / 2) : end - floor(m / 2));
        moment_out = sum(sum(k1' * k2 .* image_data));
    else
        moment_out = 0;
    end
end