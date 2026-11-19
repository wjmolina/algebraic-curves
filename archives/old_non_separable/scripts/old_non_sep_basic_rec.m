clear;
clc;

disp("getting polynomial");
n = 4;
k = 8;
[c1, a, b, c, d] = get_random_bounded_polynomial(n, 0, .5, 0, .5, 2 * k);
image_data = realize_polynomial(c1, a, b, c, d, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(image_data);

disp("down-sampling");
p = 0;
d_array = convolve_and_downsample(image_data, 2 * n - 1, 2 ^ p, 2 ^ p);
figure;
imshow(d_array / 2 ^ (2 * p));

disp("reconstructing");
c2 = non_sep_alg(d_array, n);
reconstruction = non_sep_bern_img(c2, n, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(reconstruction);

figure;
imshow(abs(image_data - reconstruction));