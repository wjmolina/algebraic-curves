function c_array = proto_no_kernel_reconstruction(image_data, n)
    m_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
    row = 1;
    for r = 0 : ceil(n / 2)
        for s = 0 : ceil(n / 2)
            col = 1;
            for k = 0 : n
                for l = 0 : n - k
                    m_matrix(row + 0, col) = (k + r) * proto_no_kernel_moment(image_data, k + r - 1, l + s);
                    m_matrix(row + 1, col) = (l + s) * proto_no_kernel_moment(image_data, k + r, l + s - 1);
                    col = col + 1;
                end
            end
            row = row + 2;
        end
    end
    % With the first coefficient fixed at one, solve for the remaining coefficients.
    c_array = [1; m_matrix(:, 2:end) \ (-m_matrix(:, 1))].';
end