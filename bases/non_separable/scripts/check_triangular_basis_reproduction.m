% Reproduce each triangular Bernstein polynomial from the same spline kernel.
n = 4;
domain_length = 12;
coefficients = non_separable_expansion_coefficients(n, domain_length);
coefficients = coefficients(:, :, :, domain_length + 1 : end, domain_length + 1 : end);
[x_grid, y_grid] = meshgrid(0 : domain_length);
kernel = wjm_b_spline(n, 1);

for k = 0 : n
    for i = 0 : k
        for j = 0 : k - i
            reconstructed = conv2(kernel, kernel, squeeze(coefficients(k + 1, i + 1, j + 1, :, :)), 'valid');
            expected = nchoosek(k, i) * nchoosek(k - i, j) .* (x_grid / domain_length) .^ i .* (y_grid / domain_length) .^ j .* (1 - (x_grid + y_grid) / domain_length) .^ (k - i - j);
            difference = abs(reconstructed - expected);
            assert(max(difference(:)) < 1e-10, 'Reproduction failed for degree %d and indices (%d, %d).', k, i, j);
        end
    end
end
