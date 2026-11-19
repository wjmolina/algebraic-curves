% clear;
% clc;

% image_data = get_image_of_power(get_random_power(4, - 1, 1, - 1, 1, false), - 1, 1, - 1, 1, 2 ^ 9 + 1, 2 ^ 9 + 1);

p_array = [6435/8388608, -7293/8388608, -8415/1048576, 9945/1048576, 85085/2097152, -109395/2097152, -153153/1048576, 255255/1048576, 3828825/4194304, 3828825/4194304, 255255/1048576, -153153/1048576, -109395/2097152, 85085/2097152, 9945/1048576, -8415/1048576, -7293/8388608, 6435/8388608];
n = 4;
b_matrix = scaling_integers(p_array);
d_array = im2double(conv2(b_matrix, b_matrix, image_data));

[height, width] = size(image_data);
x_coefficients = zeros(n + ceil(n / 2) + 1, size(d_array, 2));
y_coefficients = zeros(n + ceil(n / 2) + 1, size(d_array, 1));
for k = 0 : n + ceil(n / 2)
    x_values = pseudosplines_pol_coeffs(n, k, p_array, -1, width);
    y_values = pseudosplines_pol_coeffs(n, k, p_array, -1, height);
    x_coefficients(k + 1, :) = x_values(1 : size(d_array, 2));
    y_coefficients(k + 1, :) = y_values(1 : size(d_array, 1));
end

m_matrix = zeros(2 * (ceil(n / 2) + 1) ^ 2, nchoosek(n + 2, 2));
row = 1;

for r = 0 : ceil(n / 2)
    for s = 0 : ceil(n / 2)
        col = 1;
        
        for i = 0 : n
            for j = 0 : n - i
                if 0 <= i + r - 1
                    m_matrix(row + 0, col) = (i + r) * sum(sum(y_coefficients(j + s + 1, :)' * x_coefficients(i + r, :) .* d_array));
                else
                    m_matrix(row + 0, col) = 0;
                end
                
                if 0 <= j + s - 1
                    m_matrix(row + 1, col) = (j + s) * sum(sum(y_coefficients(j + s, :)' * x_coefficients(i + r + 1, :) .* d_array));
                else
                    m_matrix(row + 1, col) = 0;
                end
                
                col = col + 1;
            end
        end
        
        row = row + 2;
    end
end

% With the first coefficient fixed at one, solve for the remaining coefficients.
coefficients = [1; m_matrix(:, 2:end) \ (-m_matrix(:, 1))].';
reconstruction = get_image_of_power(coefficients, 0, width, 0, height, width, height);
imshow(abs(image_data - reconstruction));