clear;
clc;

% Historical inputs and outputs are relative to the original experiment root.
experiment_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
output_dir = fullfile(experiment_root, 'statistics');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

for i = 70 : - 5 : 30
    j = randi(20);
    image_data = load_noise_image(experiment_root, j);
    sample_data = awgn(convolve_and_downsample(image_data, 2 * 4 - 1, 1, 1), i);
    reconstruction = non_sep_bern_img(non_sep_alg(sample_data, 4), 4, size(sample_data, 1), size(sample_data, 2));
    denoised_reconstruction = non_sep_bern_img(non_sep_alg(wiener2(sample_data), 4), 4, size(sample_data, 1), size(sample_data, 2));
    save(fullfile(output_dir, "rec_" + j + "_" + i + ".mat"), "image_data", "sample_data", "reconstruction", "denoised_reconstruction");
end