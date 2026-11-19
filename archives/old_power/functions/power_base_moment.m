function power_base_moment_out = power_base_moment(n, i, j, sample_data, l)
if 0 <= i && 0 <= j
    k1 = dual_expansion_coefficients(n, size(sample_data, 1), l);
    k2 = dual_expansion_coefficients(n, size(sample_data, 2), l);
    power_base_moment_out = sum(sum(k1(1 + j, 1 + floor(n / 2) : end - floor(n / 2))' * k2(1 + i, 1 + floor(n / 2) : end - floor(n / 2)) .* sample_data));
else
    power_base_moment_out = 0;
end
end