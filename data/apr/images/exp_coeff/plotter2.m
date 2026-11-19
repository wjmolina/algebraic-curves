clear;
clc;

n = 4;
domain_length = 30;
c_array = non_separable_expansion_coefficients(n, domain_length);
c_c = (size(c_array, 4) + 1) / 2;

for i = 0 : n
    for j = 0 : n - i
        figure;
        image_data = reshape(c_array(1 + n, 1 + i, 1 + j, :, :), size(c_array, 4), size(c_array, 5));
        imagesc(image_data(c_c : c_c + 15, c_c : c_c + 15));
        
        xticks(1 : 16);
        xticklabels(0 : 15);
        yticks(1 : 16);
        yticklabels(0 : 15);
        
        set(gca,'YDir','normal');
        colorbar;
        title('Non-Separable Bernstein Expansion Coefficients \eta_{k, l}^{4, ' + string(i) + ', ' + string(j) + ', 4, 4, 30}');
        xlabel('k');
        ylabel('l');
    end
end