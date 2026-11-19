function image_data = get_image_of_non_separable(c_array, domain_length, x_n, y_n, c)
    [x_grid, y_grid] = meshgrid(linspace(0, domain_length, x_n), linspace(0, domain_length, y_n));
    image_data = zeros(size(x_grid));
    n = (sqrt(8 * numel(c_array) + 1) - 3) / 2;
    if ~isvector(c_array) || n < 0 || n ~= floor(n)
        error('Expected a coefficient vector of length (n + 1) * (n + 2) / 2.');
    end
    k = 1;
    new_l = c * domain_length;
    for i = 0 : n
        for j = 0 : n - i
            image_data = image_data + c_array(k) * new_l ^ (- n) * nchoosek(n, i) * nchoosek(n - i, j) * x_grid .^ i .* y_grid .^ j .* (new_l - x_grid - y_grid) .^ (n - i - j);
            k = k + 1;
        end
    end
    image_data = im2double(image_data <= 0);
end