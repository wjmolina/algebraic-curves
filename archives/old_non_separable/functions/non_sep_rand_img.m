function [non_sep_rand_img_out, coefficient_array] = non_sep_rand_img(n, a, b)
    validateattributes(n, {'numeric'}, {'scalar', 'real', 'finite', 'integer', '>=', 2});
non_sep_rand_img_out = zeros(a, b);
while non_sep_rand_img_check(non_sep_rand_img_out)
    non_sep_rand_img_out = zeros(a, b);
    coefficient_array = zeros(1, nchoosek(n + 2, 2));
    k = 1;
    for i = 0 : n
        for j = 0 : n - i
            coeff = (2 * rand() - 1) * 100;
            coefficient_array(k) = coeff;
            k = k + 1;
            non_sep_rand_img_out = non_sep_rand_img_out + coeff * non_sep_bern(1, a, b, n, i, j);
        end
    end
    non_sep_rand_img_out = im2double(non_sep_rand_img_out <= 0);
end
end