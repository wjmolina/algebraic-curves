function non_sep_moment_out = non_sep_moment(image_data, n, i, j)
if 0 <= n && 0 <= i && 0 <= j && i + j <= n
    k_array = non_sep_exp_coeffs(n, size(image_data, 1), size(image_data, 2));
    non_sep_moment_out = sum(sum(reshape(k_array(1 + i, 1 + j, 1 + floor(n / 2) : end - floor(n / 2), 1 + floor(n / 2) : end - floor(n / 2)), size(image_data, 1), size(image_data, 2)) .* image_data));
else
    non_sep_moment_out = 0;
end
end