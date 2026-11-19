% Display the reproduction error for the degree-four basis functions.
clear;
clc;

n = 4;
domain_length = 12;
coefficients = non_separable_expansion_coefficients(n, domain_length);
coefficients = coefficients(:, :, :, domain_length + 1 : end, domain_length + 1 : end);
[x_grid, y_grid] = meshgrid(0 : domain_length);
kernel = wjm_b_spline(n, 1);
for i = 0 : n
    for j = 0 : n - i
        reconstructed = conv2(kernel, kernel, squeeze(coefficients(n + 1, i + 1, j + 1, :, :)), 'valid');
        expected = nchoosek(n, i) * nchoosek(n - i, j) .* (x_grid / domain_length) .^ i .* (y_grid / domain_length) .^ j .* (1 - (x_grid + y_grid) / domain_length) .^ (n - i - j);
        surf(abs(reconstructed - expected));
        pause(.5);
    end
end
