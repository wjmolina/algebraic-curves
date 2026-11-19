function [c_array, d_array] = proto_power_reconstruction_db(image_data, n, target_snr_db, t)
    [k_array, s_array] = get_daub_valid_conv_coeff(n + ceil(n / 2) + 1, 1 - 100, size(image_data, 1) + 100);
%     s_array = scaling_integers(2 * dbwavf("db" + (n + ceil(n / 2) + 1)));
    d_array = conv2(s_array, s_array, image_data);
    if target_snr_db ~= - 1
        d_array = imnoise(d_array, 'gaussian', 0, 10 ^ (- target_snr_db / 10) * mean(mean(d_array .^ 2)));
    end
%     k_array = proto_power_expansion_coefficients_db(n + ceil(n / 2) + 1, - length(d_array) - 50, length(d_array) + 50);
    m_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;
    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;
            for i = 0 : n
                for j = 0 : n - i
                    if 0 <= i + r - 1
                        m_matrix(row + 0, col) = (i + r) * sum(sum(k_array(j + s + 1, t : t + length(d_array) - 1)' * k_array(i + r + 0, t : t + length(d_array) - 1) .* d_array));
                    else
                        m_matrix(row + 0, col) = 0;
                    end
                    if 0 <= j + s - 1
                        m_matrix(row + 1, col) = (j + s) * sum(sum(k_array(j + s + 0, t : t + length(d_array) - 1)' * k_array(i + r + 1, t : t + length(d_array) - 1) .* d_array));
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