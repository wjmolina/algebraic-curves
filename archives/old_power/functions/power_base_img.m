function power_base_img_out = power_base_img(k_array, n, r, s, l)
image_data = zeros(r, s);
k = 1;
for i = 0 : n
    for j = 0 : n - i
        image_data = image_data + k_array(k) .* power_base(i, j, r, s, l);
        k = k + 1;
    end
end
power_base_img_out = im2double(image_data <= 0);
end