function power_base_alg_out = power_base_alg(sample_data, n, l)
a_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
row = 1;
for r = 0 : ceil(n / 2)
    for s = 0 : ceil(n / 2)
        col = 1;
        for i = 0 : n
            for j = 0 : n - i
                a_matrix(row + 0, col) = (i + r) * power_base_moment(n + ceil(n / 2), i + r - 1, j + s, sample_data, l);
                a_matrix(row + 1, col) = (j + s) * power_base_moment(n + ceil(n / 2), i + r, j + s - 1, sample_data, l);
                col = col + 1;
            end
        end
        row = row + 2;
    end
end
power_base_alg_out = transpose(lsqlin(a_matrix, zeros(size(a_matrix, 1), 1), [], [], eye(1, size(a_matrix, 2)), 1, [], []));
% power_base_alg_out = transpose(quadprog(transpose(a_matrix) * a_matrix, [], [], [], eye(1, size(a_matrix, 2)), 1));
end