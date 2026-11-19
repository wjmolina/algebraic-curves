%          n = max degree of poly
% a, b, c, d = bounds of domain
%        pad = clearance of image with border

function [coeffs, a, b, c, d] = get_random_bounded_polynomial(n, a, b, c, d, pad)
    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'integer', '>=', 2});
    [x_grid, y_grid] = meshgrid(linspace(a, b, 2 ^ 9), linspace(c, d, 2 ^ 9));
    img = zeros(size(x_grid));
    while img(1, 1) == 1 || length(unique(img(1 : pad, :))) ~= 1 || length(unique(img(end + 1 - pad : end, :))) ~= 1 || length(unique(img(:, 1 : pad))) ~= 1 || length(unique(img(:, end + 1 - pad : end))) ~= 1 || length(unique(img)) ~= 2
        img = zeros(size(x_grid));
        coeffs = zeros(1, nchoosek(n + 2, 2));
        cnt = 1;
        for i = 0 : n
            for j = 0 : n - i
                coeffs(cnt) = 2 * rand() - 1;
                img = img + coeffs(cnt) * x_grid .^ i .* y_grid .^ j;
                cnt = cnt + 1;
            end
        end
        img = im2double(img <= 0);
    end
end