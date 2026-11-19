clear;
clc;

n = 4;
a = 128;
b = 128;

nsec = non_sep_exp_coeffs(n, a, b);

for i = 0 : n
    for j = 0 : n - i
        imshow(conv2(wjm_b_spline(n, 1), wjm_b_spline(n, 1), reshape(nsec(1 + i, 1 + j, :, :), size(nsec, 3), size(nsec, 4)), 'valid'));
        pause(.25);
    end
end