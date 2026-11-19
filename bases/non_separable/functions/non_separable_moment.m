function m_matrix = non_separable_moment(n, i, j, sample_data, c_array, c_x, c_y)
    if n >= 0 && i >= 0 && j >= 0 && i + j <= n
        m_matrix = sum(sum(reshape(c_array(1 + n, 1 + i, 1 + j, c_y, c_x), size(c_y, 2), size(c_x, 2)) .* sample_data));
    else
        m_matrix = 0;
    end
end