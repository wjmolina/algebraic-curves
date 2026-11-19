function image_data = get_image_of_power(c_array, x_a, x_b, y_a, y_b, x_n, y_n)
    [x_grid, y_grid] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    image_data = zeros(size(x_grid));
    n = (sqrt(8 * numel(c_array) + 1) - 3) / 2;
    if ~isvector(c_array) || n < 0 || n ~= floor(n)
        error('Expected a coefficient vector of length (n + 1) * (n + 2) / 2.');
    end
    k = 1;
    for i = 0 : n
        for j = 0 : n - i
            image_data = image_data + c_array(k) * x_grid .^ i .* y_grid .^ j;
            k = k + 1;
        end
    end
    image_data = im2double(image_data <= 0);
end