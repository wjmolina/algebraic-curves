function m_matrix = get_random_symmetric_positive_definite_matrix()
    v_matrix = zeros(3, 3);

    v_matrix(1, :) = (2 * rand() - 1) * [1, 1, 1];
    v_matrix(2, :) = (2 * rand() - 1) * [1, 1, 0];
    v_matrix(3, :) = (2 * rand() - 1) * [1, 0, 0];

    m_matrix = zeros(3, 3);

    for i = 1 : 3
        for j = 1 : 3
            m_matrix(i, j) = dot(v_matrix(i, :), v_matrix(j, :));
        end
    end
end