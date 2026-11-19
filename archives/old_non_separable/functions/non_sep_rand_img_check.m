function non_sep_rand_img_check_out = non_sep_rand_img_check(image_data)
non_sep_rand_img_check_out = false;
if length(unique(union(union(union(image_data(1, :), image_data(end, :)), image_data(:, 1)), image_data(:, end)))) ~= 1 || length(unique(image_data)) == 1 || length(unique(triu(rot90(image_data)))) ~= 1
    non_sep_rand_img_check_out = true;
end
end