function [p_array, image_data] = get_random_power(n, x_a, x_b, y_a, y_b, is_other_mode)
    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'integer', '>=', 2});
    [x_grid, y_grid] = meshgrid(linspace(x_a, x_b, 2 ^ 9 + 1), linspace(y_a, y_b, 2 ^ 9 + 1));
    image_data = zeros(size(x_grid));
    while image_data(1, 1) == 1 || length(unique(image_data(1, :))) ~= 1 || length(unique(image_data(end, :))) ~= 1 || length(unique(image_data(:, 1))) ~= 1 || length(unique(image_data(:, end))) ~= 1 || length(unique(image_data)) ~= 2
        image_data = zeros(size(x_grid));
        p_array = zeros(1, nchoosek(n + 2, 2));
        k = 1;
        if ~ is_other_mode
            for i = 0 : n
                for j = 0 : n - i
                    p_array(k) = 2 * rand() - 1;
                    image_data = image_data + p_array(k) * x_grid .^ i .* y_grid .^ j;
                    k = k + 1;
                end
            end
        else
            for l = 0 : n
                for i = 0 : l
                    j = l - i;
                    p_array(k) = 2 * rand() - 1;
                    image_data = image_data + p_array(k) * x_grid .^ i .* y_grid .^ j;
                    k = k + 1;
                end
            end
        end
        image_data = im2double(image_data <= 0);
    end
end