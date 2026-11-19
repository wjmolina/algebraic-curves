function bernstein_factor_out = separable_factor(n, i, r)
    bernstein_factor_out = 1;
    if r > 0
        for k = 1 : r
            bernstein_factor_out = bernstein_factor_out * (i + k) / (n + k);
        end
    end
end