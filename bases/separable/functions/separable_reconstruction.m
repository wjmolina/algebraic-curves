function sep_algorithm_out = separable_reconstruction(image_data, n)
    % System of separable_moment Equations
    a_matrix = zeros(2 * (n + 1) ^ 2, (n + 1) ^ 2);
    row = 1;
    for r = 0 : n
        for s = 0 : n
            for i = 0 : n
                for j = 0 : n
                    helper = separable_factor(n, i, r) * separable_factor(n, j, s);
                    a_matrix(row + 0, 1 + (n + 1) * i + j) = helper * (n + r) * (separable_moment(image_data, n + r - 1, n + s, i + r - 1, j + s) - separable_moment(image_data, n + r - 1, n + s, i + r, j + s));
                    a_matrix(row + 1, 1 + (n + 1) * i + j) = helper * (n + s) * (separable_moment(image_data, n + r, n + s - 1, i + r, j + s - 1) - separable_moment(image_data, n + r, n + s - 1, i + r, j + s));
                end
            end
            row = row + 2;
        end
    end
    % Power Basis Constraints
    b_matrix = zeros(1, size(a_matrix, 2));
    b_matrix(1, 1) = 1;
    row = 2;
    for r = 0 : n
        for s = 0 : n
            if r + s > n
                b_matrix(row, :) = zeros(1, size(a_matrix, 2));
                for i = 0 : r
                    for j = 0 : s
                        b_matrix(row, 1 + (n + 1) * i + j) = nchoosek(n, r) * nchoosek(n, s) * nchoosek(r, i) * nchoosek(s, j) * (- 1) ^ (r + s - i - j);
                    end
                end
                row = row + 1;
            end
        end
    end
    sep_algorithm_out = lsqlin(a_matrix, zeros(size(a_matrix, 1), 1), [], [], b_matrix, eye(size(b_matrix, 1), 1), [], [], []);
end