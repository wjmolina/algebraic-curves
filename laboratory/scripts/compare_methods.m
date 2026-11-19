clear;
clc;

n = 4;
domain_length = 20;
m_x = 25;
m_y = 25;
target_snr_db = 50;
n_total = 100;

x_n = m_x * domain_length + 1;
y_n = m_y * domain_length + 1;
d_array = zeros(n_total, 4);

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
load(fullfile(repo_root, 'data', 'statistics', 'images_02.mat'), 'image_batch');
triangle_scale = 2;

for i = 1 : n_total
    image_data = reshape(image_batch(i, :, :), [size(image_batch, 2), size(image_batch, 3)]);
    
    d_array(i, 1) = sorensen_dice_coefficient(image_data, get_image_of_power(power_reconstruction(image_data, n, 0, domain_length, 0, domain_length, target_snr_db), 0, domain_length, 0, domain_length, x_n, y_n));
    
    d_array(i, 2) = sorensen_dice_coefficient(image_data, get_image_of_non_separable(non_separable_reconstruction(image_data, n, domain_length, target_snr_db, triangle_scale), domain_length, x_n, y_n, triangle_scale));
    
    [qp_e, qp] = vetterli_reconstruction(image_data, 0, domain_length, 0, domain_length, x_n, y_n, target_snr_db);
    d_array(i, 3) = sorensen_dice_coefficient(image_data, qp);
    d_array(i, 4) = sorensen_dice_coefficient(image_data, qp_e);
    
    fprintf('%d/%d\n', i, n_total);
end

% figure;
% histogram(s_array(:, k));

figure;
boxplot(d_array, 'Labels', {'Power', 'Triangular Bernstein', 'Vetterli QP', 'Vetterli enhanced'}, 'Colors', 'k', 'Whisker', 1000000000);
set(gca, 'FontSize', 9);
title(sprintf('Reconstruction comparison of %d images at %g dB', n_total, target_snr_db));
ylabel('Sørensen–Dice Coefficient');