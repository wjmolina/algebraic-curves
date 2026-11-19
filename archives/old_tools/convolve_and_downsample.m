function convolve_and_downsample_out = convolve_and_downsample(image_data, n, m, t)
    helper = conv2(wjm_b_spline(n, m), wjm_b_spline(n, m), image_data, 'same');
    convolve_and_downsample_out = helper(1 : t : end, 1 : t : end);
end