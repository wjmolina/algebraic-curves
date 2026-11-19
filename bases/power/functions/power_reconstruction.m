function [c_array, d_array] = power_reconstruction(image_data, n, x_a, x_b, y_a, y_b, target_snr_db)
    m_x = (size(image_data, 2) - 1) / (x_b - x_a);
    m_y = (size(image_data, 1) - 1) / (y_b - y_a);

    if floor(m_y) ~= m_y || floor(m_x) ~= m_x
        disp("image domain and resolution inconsistency"); % This can be fixed by up-sampling and shifting splines, I think.
        return;
    end
    
    d_t = conv2(wjm_b_spline(n + ceil(n / 2), m_y), wjm_b_spline(n + ceil(n / 2), m_x), image_data) / (m_x * m_y);
    d_t_rx = (size(d_t, 2) - size(image_data, 2)) / 2 + 1;
    d_t_ry = (size(d_t, 1) - size(image_data, 1)) / 2 + 1;

    d_array = d_t(d_t_ry - floor((d_t_ry - 1) / m_y) * m_y : m_y : end, d_t_rx - floor((d_t_rx - 1) / m_x) * m_x : m_x : end);
    d_rx = (size(d_array, 2) - size(image_data(1 : m_y : end, 1 : m_x : end), 2)) / 2 + 1;
    d_ry = (size(d_array, 1) - size(image_data(1 : m_y : end, 1 : m_x : end), 1)) / 2 + 1;
    
    if target_snr_db ~= - 1
        d_array = imnoise(d_array, 'gaussian', 0, 10 ^ (- target_snr_db / 10) * mean(mean(d_array .^ 2)));
    end

    l = max(abs([1 - d_rx + x_a, 1 - d_ry + y_a, size(d_array, 2) - d_rx + x_a, size(d_array, 1) - d_ry + y_a]));
    
    k_array = power_expansion_coefficients(n + ceil(n / 2), - l, l);
    k_r = (size(k_array, 2) + 1) / 2;
    k_x = k_r + 1 - d_rx + x_a : k_r + size(d_array, 2) - d_rx + x_a;
    k_y = k_r + 1 - d_ry + y_a : k_r + size(d_array, 1) - d_ry + y_a;

    m_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;

    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;

            for i = 0 : n
                for j = 0 : n - i
                    if 0 <= i + r - 1
                        m_matrix(row + 0, col) = (i + r) * sum(sum(k_array(j + s + 1, k_y)' * k_array(i + r + 0, k_x) .* d_array));
                    else
                        m_matrix(row + 0, col) = 0;
                    end

                    if 0 <= j + s - 1
                        m_matrix(row + 1, col) = (j + s) * sum(sum(k_array(j + s + 0, k_y)' * k_array(i + r + 1, k_x) .* d_array));
                    else
                        m_matrix(row + 1, col) = 0;
                    end

                    col = col + 1;
                end
            end

            row = row + 2;
        end
    end

    % With the first coefficient fixed at one, solve for the remaining coefficients.
    c_array = [1; m_matrix(:, 2:end) \ (-m_matrix(:, 1))].';
end