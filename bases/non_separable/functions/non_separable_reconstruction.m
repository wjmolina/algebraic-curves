function [c_array, d_array] = non_separable_reconstruction(image_data, n, domain_length, target_snr_db, c)
    m_x = (size(image_data, 2) - 1) / domain_length;
    m_y = (size(image_data, 1) - 1) / domain_length;

    if floor(m_y) ~= m_y || floor(m_x) ~= m_x
        disp("image domain and resolution inconsistency");
        return;
    end
    
    kernel_degree = n + 2 * ceil(n / 2) - 1;
    d_t = conv2(wjm_b_spline(kernel_degree, m_y), wjm_b_spline(kernel_degree, m_x), image_data) / (m_x * m_y);
    d_t_x = (size(d_t, 2) - size(image_data, 2)) / 2 + 1;
    d_t_y = (size(d_t, 1) - size(image_data, 1)) / 2 + 1;

    d_array = d_t(d_t_y - floor((d_t_y - 1) / m_y) * m_y : m_y : end, d_t_x - floor((d_t_x - 1) / m_x) * m_x : m_x : end);
    d_x = (size(d_array, 2) - size(image_data(1 : m_y : end, 1 : m_x : end), 2)) / 2 + 1;
    d_y = (size(d_array, 1) - size(image_data(1 : m_y : end, 1 : m_x : end), 1)) / 2 + 1;
    
    if target_snr_db ~= - 1
        d_array = imnoise(d_array, 'gaussian', 0, 10 ^ (- target_snr_db / 10) * mean(mean(d_array .^ 2)));
    end
    
    new_l = c * domain_length;
    
    k_array = non_separable_expansion_coefficients(kernel_degree, new_l);
    k_c = (size(k_array, 4) + 1) / 2;
    k_x = k_c + 1 - d_x : k_c + size(d_array, 2) - d_x;
    k_y = k_c + 1 - d_y : k_c + size(d_array, 1) - d_y;

    a_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for a = 0 : ceil(n / 2)
        for b = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    a_matrix(row + 0, col) = (non_separable_moment(n + a + b - 1, i + a - 1, j + b + 0, d_array, k_array, k_x, k_y) - non_separable_moment(n + a + b - 1, i + a, j + b, d_array, k_array, k_x, k_y)) * n * non_separable_factor(n, a, b, i, j, new_l) / new_l;
                    a_matrix(row + 1, col) = (non_separable_moment(n + a + b - 1, i + a + 0, j + b - 1, d_array, k_array, k_x, k_y) - non_separable_moment(n + a + b - 1, i + a, j + b, d_array, k_array, k_x, k_y)) * n * non_separable_factor(n, a, b, i, j, new_l) / new_l;

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end
    
    % With the first coefficient fixed at one, solve for the remaining coefficients.
    c_array = [1; a_matrix(:, 2:end) \ (-a_matrix(:, 1))].';
end