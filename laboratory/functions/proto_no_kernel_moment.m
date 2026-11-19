function m_matrix = proto_no_kernel_moment(image_data, i, j)
    if i >= 0 && j >= 0
        [x_grid, y_grid] = meshgrid(0 : size(image_data, 2) - 1, 0 : size(image_data, 1) - 1);
        m_matrix = sum(sum(x_grid .^ i .* y_grid .^ j .* image_data));
    else
        m_matrix = 0;
    end
end