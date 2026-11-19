function image_data = get_random_power_stably_bounded(x_a, x_b, y_a, y_b, x_n, y_n)
    a_matrix = get_random_symmetric_positive_definite_matrix();
    [x_grid, y_grid] = meshgrid(linspace(x_a, x_b, x_n), linspace(y_a, y_b, y_n));
    while true
        c_array = 2 * rand(1, 15) - 1;
        c_array(1) = 1;
        image_data = arrayfun(@(x, y) [x ^ 2, x * y, y ^ 2] * a_matrix * [x ^ 2, x * y, y ^ 2]', x_grid, y_grid);
        cnt = 1;
        for i = 0 : 4
            for j = 0 : 4 - i
                if i + j < 4
                    image_data = image_data + c_array(cnt) * x_grid .^ i .* y_grid .^ j;
                else
                    c_array(cnt) = 0;
                end
                cnt = cnt + 1;
            end
        end
        image_data = im2double(image_data <= 0);
        if ~ (image_data(1, 1) == 1 || length(unique(image_data(1, :))) ~= 1 || length(unique(image_data(end, :))) ~= 1 || length(unique(image_data(:, 1))) ~= 1 || length(unique(image_data(:, end))) ~= 1 || length(unique(image_data)) ~= 2)
            break;
        end
    end
end