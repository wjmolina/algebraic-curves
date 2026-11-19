clear;
clc;

stats = zeros(100, 50);

for image_n = 1 : 100
    disp(image_n);
    [~, image_data] = get_random_non_separable(4, 3);
    
    for triangle_scale = 1 : 50
        reconstruction = get_image_of_non_separable(non_separable_reconstruction(image_data, 4, 1, - 1, triangle_scale), 1, 257, 257, triangle_scale);
        stats(image_n, triangle_scale) = sorensen_dice_coefficient(image_data, reconstruction);
    end
end

boxplot(stats, 'whisker', 1000);
xlabel('Triangle scale c');
ylabel('Sørensen–Dice Coefficient');
title('Bernstein Reconstruction of 100 Algebraic Shapes in [0,1]^2');