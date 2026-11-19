clear;
clc;

experiment_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
output_dir = fullfile(experiment_root, 'statistics');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

for j = 25 : - 5 : 25
    error_n = [];
    error_d = [];
    for i = 1 : 20
        image_data = load_noise_image(experiment_root, i);
        sample_data = awgn(convolve_and_downsample(image_data, 2 * 4 - 1, 1, 1), j);
        k_array = abs(image_data - non_sep_bern_img(non_sep_alg(sample_data, 4), 4, size(image_data, 1), size(image_data, 2)));
        error_n = [error_n, sum(sum(k_array)) / numel(k_array)];
        k_array = abs(image_data - non_sep_bern_img(non_sep_alg(wiener2(sample_data), 4), 4, size(image_data, 1), size(image_data, 2)));
        error_d = [error_d, sum(sum(k_array)) / numel(k_array)];
    end
    save(fullfile(output_dir, "error_" + j + ".mat"), "error_n", "error_d");
end