function [p_array, image_data] = get_random_separable(n, x_a, x_b, y_a, y_b)
    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'integer', '>=', 2});
    [x_grid, y_grid] = meshgrid(linspace(x_a, x_b, 2 ^ 9 + 1), linspace(y_a, y_b, 2 ^ 9 + 1));
    image_data = zeros(size(x_grid));
    while image_data(1, 1) == 1 || length(unique(image_data(1, :))) ~= 1 || length(unique(image_data(end, :))) ~= 1 || length(unique(image_data(:, 1))) ~= 1 || length(unique(image_data(:, end))) ~= 1 || length(unique(image_data)) ~= 2
        image_data = zeros(size(x_grid));
        p_array = zeros(1, (n + 1) ^ 2);
        k = 1;
        for i = 0 : n
            for j = 0 : n
                p_array(k) = 2 * rand() - 1;
                image_data = image_data + p_array(k) .* nchoosek(n, i) .* ((x_grid - x_a) / (x_b - x_a)) .^ i .* (1 - (x_grid - x_a) / (x_b - x_a)) .^ (n - i) .* nchoosek(n, j) .* ((y_grid - y_a) / (y_b - y_a)) .^ j .* (1 - (y_grid - y_a) / (y_b - y_a)) .^ (n - j);
                k = k + 1;
            end
        end
        image_data = im2double(image_data <= 0);
    end
end