function display_separable_to_power_coefficients(bernstein_coefficients)
    width = sqrt(numel(bernstein_coefficients));
    if width < 1 || width ~= floor(width)
        error('Expected a nonzero square number of Bernstein coefficients.');
    end
    n = width - 1;
    power_coefficients = zeros(1, numel(bernstein_coefficients));
    index = 1;
    for r = 0 : n
        for s = 0 : n
            helper = 0;
            for i = 0 : r
                for j = 0 : s
                    helper = helper + bernstein_coefficients(1 + width * i + j) * nchoosek(n, r) * nchoosek(n, s) * nchoosek(r, i) * nchoosek(s, j) * (- 1) ^ (r + s - i - j);
                end
            end
            power_coefficients(index) = helper;
            index = index + 1;
        end
    end
    for i = 0 : n
        for j = 0 : n
            disp([num2str(power_coefficients(1 + width * i + j)), ' * x ^ ', num2str(i), ' * y ^ ', num2str(j)]);
        end
    end
end