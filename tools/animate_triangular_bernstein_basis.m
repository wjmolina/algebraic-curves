% Animate the degree-n Bernstein basis functions on the closed triangle.

clear;
clc;

n = 4;
domain_length = 2;
k_array = linspace(0, domain_length, 512);

[xx, yy] = meshgrid(k_array, k_array);

k = 1;

for i = 0 : n
    for j = 0 : n - i
        f = figure;        
                image_data = nchoosek(n, i) .* nchoosek(n - i, j) .* xx .^ i .* yy .^ j .* (domain_length - xx - yy) .^ (n - i - j) ./ domain_length .^ n;
                imagesc(k_array, k_array, image_data .* (xx + yy <= domain_length));
                title(sprintf("B_{%i,%i}^{%i} Restricted to Triangle", i, j, n));        
%         image_data = nchoosek(n, i) .* nchoosek(n - i, j) .* xx .^ i .* yy .^ j .* (2 * domain_length - xx - yy) .^ (n - i - j) ./ (2 * domain_length) .^ n;
%         imagesc(k_array, k_array, image_data);
%         title(sprintf("B_{%i,%i}^{%i} Restricted to Square", i, j, n));
        
        ax = gca;
        ax.YDir = "normal";
        colorbar;
        caxis([0, 1]);
        xticks(0 : .2 : 2);
        yticks(0 : .2 : 2);
        
        [x_grid, map] = rgb2ind(frame2im(getframe(f)), 256);        
        if k == 1
            imwrite(x_grid, map, 'non_sep_trg.gif', 'gif', 'Loopcount', inf);
        else
            imwrite(x_grid, map, 'non_sep_trg.gif', 'gif', 'WriteMode', 'append');
        end        
        k = k + 1;
    end
end