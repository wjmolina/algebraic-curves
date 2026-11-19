function [p_array, image_data] = get_random_non_separable(n, domain_length)
    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'integer', '>=', 2});
    [x_grid, y_grid] = meshgrid(linspace(0, domain_length, 2 ^ 8 + 1), linspace(0, domain_length, 2 ^ 8 + 1));
    image_data = zeros(size(x_grid));
    while image_data(1, 1) == 1 || length(unique(image_data(1, :))) ~= 1 || length(unique(image_data(end, :))) ~= 1 || length(unique(image_data(:, 1))) ~= 1 || length(unique(image_data(:, end))) ~= 1 || length(unique(image_data)) ~= 2
        image_data = zeros(size(x_grid));
        p_array = zeros(1, nchoosek(n + 2, 2));
        k = 1;
        for i = 0 : n
            for j = 0 : n - i
                p_array(k) = 2 * rand() - 1;
                image_data = image_data + p_array(k) * nchoosek(n, i) * nchoosek(n - i, j) * x_grid .^ i .* y_grid .^ j .* (2 * domain_length - x_grid - y_grid) .^ (n - i - j) / (2 * domain_length) ^ n;
                k = k + 1;
            end
        end
        image_data = im2double(image_data <= 0);
    end
end