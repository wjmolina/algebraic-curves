function get_separable_bernstein_image_out = non_sep_bern_img(k_array, n, a, b)
image_data = zeros(a, b);
k = 1;
for i = 0 : n
    for j = 0 : n - i
        image_data = image_data + k_array(k) .* non_sep_bern(1, a, b, n, i, j);
        k = k + 1;
    end
end
get_separable_bernstein_image_out = im2double(image_data <= 0);
end