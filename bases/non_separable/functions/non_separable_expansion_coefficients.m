function n_total = non_separable_expansion_coefficients(n, domain_length)
    p_array = power_expansion_coefficients(n, - domain_length, domain_length);
    n_total = zeros(n + 1, n + 1, n + 1, size(p_array, 2), size(p_array, 2));
    
    for k = 0 : n
        for i = 0 : k
            for j = 0 : k - i
                s_array = zeros(size(p_array, 2), size(p_array, 2));

                for a = 0 : k - i - j
                    for b = 0 : k - i - j - a
                        c = k - i - j - a - b;
                        
                        s_array = s_array + nchoosek(k, i) * nchoosek(k - i, j) * factorial(k - i - j) / (factorial(a) * factorial(b) * factorial(c)) * domain_length ^ (a - k) * (- 1) ^ (b + c) * p_array(1 + c + j, :)' * p_array(1 + b + i, :);
                    end
                end

                n_total(1 + k, 1 + i, 1 + j, :, :) = s_array;
            end
        end
    end
end