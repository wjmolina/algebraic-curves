function non_sep_algorithm_out = non_sep_alg(image_data, n)
a_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
row = 1;
for r = 0 : ceil(n / 2)
    for s = 0 : ceil(n / 2)
        col = 1;
        for i = 0 : n
            for j = 0 : n - i
                helper_1 = n * non_sep_bern_prod(n, i, j, r, s);
                helper_2 = non_sep_moment(image_data, n + r + s - 1, i + r, j + s);
                a_matrix(row + 0, col) = helper_1 * (non_sep_moment(image_data, n + r + s - 1, i + r - 1, j + s) - helper_2);
                a_matrix(row + 1, col) = helper_1 * (non_sep_moment(image_data, n + r + s - 1, i + r, j + s - 1) - helper_2);
                col = col + 1;
            end
        end
        row = row + 2;
    end
end
% non_sep_algorithm_out = transpose(lsqlin(a_matrix, zeros(size(a_matrix, 1), 1), [], [], eye(1, size(a_matrix, 2)), 1, [], []));
% With the first coefficient fixed at one, solve for the remaining coefficients.
non_sep_algorithm_out = [1; a_matrix(:, 2:end) \ (-a_matrix(:, 1))].';
end