function wjm_b_spline_out = wjm_b_spline(n, m)
    x_grid = - ceil((n * m + m) / 2) : ceil((n * m + m) / 2);
    b1(1, :) = (- m / 2 < x_grid & x_grid < m / 2) * 1 + (- m / 2 == x_grid | x_grid == m / 2) * .5;
    b2(1, :) = (- m < x_grid & x_grid < 0) * 1 + (- m == x_grid | x_grid == 0) * .5;
    for i = 1 : n
        b1(i + 1, :) = ((x_grid / m + (i + 1) / 2) .* b2(i, :) + ((i + 1) / 2 - x_grid / m) .* [zeros(1, m), b2(i, 1 : end - m)]) / i;
        b2(i + 1, :) = ((x_grid / m + (i + 2) / 2) .* [b1(i, m + 1 : end), zeros(1, m)] + (i / 2 - x_grid / m) .* b1(i, :)) / i;
    end
    wjm_b_spline_out = b1(end, :);
    wjm_b_spline_out = wjm_b_spline_out(2 : end - 1);
end