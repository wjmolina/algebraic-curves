% INPUT
%     coeffs = coefficients of poly
% a, b, c, d = bounds of domain
%       r, s = resolution of domain
% OUTPUT
% img = r by s matrix taking values in {0, 1} of poly p^{-1}((-inf, 0]) \cap [a, b] × [c, d] induced by coeffs

function img = realize_polynomial(c_array, x_bound_a, x_bound_b, y_bound_a, y_bound_b, x_count, y_count)
    [x_grid, y_grid] = meshgrid(linspace(x_bound_a, x_bound_b, x_count), linspace(y_bound_a, y_bound_b, y_count));
    img = zeros(size(x_grid));
    cnt = 1;
    n = (sqrt(8 * numel(c_array) + 1) - 3) / 2;
    if ~isvector(c_array) || n < 0 || n ~= floor(n)
        error('Expected a coefficient vector of length (n + 1) * (n + 2) / 2.');
    end
    for i = 0 : n
        for j = 0 : n - i
            img = img + c_array(cnt) * x_grid .^ i .* y_grid .^ j;
            cnt = cnt + 1;
        end
    end
    img = im2double(img <= 0)';
end