clear;
clc;

disp("getting polynomial");
n = 4;
k = 9;
[c1, a, b, c, d] = get_random_bounded_polynomial(n, 0, 1, 0, 1, 5 * k);
image_data = realize_polynomial(c1, a, b, c, d, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(image_data);

disp("down-sampling");
p = 4;
d_array = convolve_and_downsample(image_data, n + ceil(n / 2), 2 ^ p, 2 ^ p);
figure;
imshow(d_array / 2 ^ (2 * p));

disp("reconstructing");
c2 = power_base_alg(d_array, n, 1);
reconstruction = realize_polynomial(c2, 0, 1, 0, 1, 2 ^ k + 1, 2 ^ k + 1);
figure;
imshow(reconstruction);

figure;
imshow(abs(image_data - reconstruction));