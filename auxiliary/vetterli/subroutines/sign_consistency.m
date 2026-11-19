% This subroutine emforces sign consistency with the available pixels; in
% simple words, it tries to find the closest fit to the null space of the
% moments that generates almost the same pixels





%%%%%%%%%%%%%%%%%%%%%
% Finding the polynomial degree vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% auxillary variables for storing the monomial powers
x_power     = 0;
y_power     = 0;

x_power_vec = zeros(1 , (sum_degree+1) * (sum_degree+2) / 2 );
y_power_vec = zeros(1 , (sum_degree+1) * (sum_degree+2) / 2 );



% finding the order of the monomials up to degree "sum_degree"
term        = 0;
while x_power + y_power <= sum_degree
    term    = term + 1;
    
    x_power_vec(term)   = x_power;
    y_power_vec(term)   = y_power;
    
    if y_power == 0
        y_power     = x_power + y_power + 1;
        x_power     = 0;
        
    else
        y_power     = y_power - 1;
        x_power     = x_power + 1;
        
    end

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial evaluation matrix over the sampling grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poly_eval_mat   = zeros( length(m_x)*length(m_y)  ,  (sum_degree+1) * (sum_degree+2) / 2 );
for m_x_ind = 1 : length(m_x)
    for m_y_ind = 1 : length(m_y)
        row_ind     = m_y_ind + (m_x_ind - 1) * length(m_y);
        poly_eval_mat(row_ind , :)  = (m_x(m_x_ind) .^ x_power_vec) .* (m_y(end-m_y_ind+1) .^ y_power_vec);
    end
end





%%%%%%%%%%%%%%%%%%%%%%
% Quadratic programming for initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh              = moment_matrix' * moment_matrix;
ff              = zeros((sum_degree+1) * (sum_degree+2) / 2   ,   1);

a_eq            = [1   ,   zeros(1 , (sum_degree+1) * (sum_degree+2) / 2 - 1)];
b_eq            = sign( pixelized_noisy(m_y==0 , m_x==0) - 0.5 );

x0              = poly_coefs_est';

a_ineq          = poly_eval_mat;
b_ineq          = ones(length(m_x)*length(m_y) , 1) * inf;
for m_x_ind = 1 : length(m_x)
    for m_y_ind = 1 : length(m_y)
        if pixelized_noisy(m_y_ind , m_x_ind) >= min(0.7 , max(max(pixelized_noisy)))
            row_ind     = m_y_ind + (m_x_ind - 1) * length(m_y);
            b_ineq(row_ind)         = 0;
            a_ineq(row_ind , :)     = -a_ineq(row_ind , :);
            
        elseif pixelized_noisy(m_y_ind , m_x_ind) <= 0.3
            row_ind     = m_y_ind + (m_x_ind - 1) * length(m_y);
            b_ineq(row_ind)         = 0;
            
        end
    end
end

options         = optimset('MaxIter' , 2e2 , 'Algorithm' , 'interior-point-convex');


% The initial minimization, the output of which is used as the initial
% point for the next iteration
[poly_coefs_qp , cost , exitflag]   = quadprog(hh*10000 , ff , a_ineq(b_ineq~=inf , :) , b_ineq(b_ineq~=inf) , a_eq , b_eq , [] , [] , x0 , options);
%[poly_coefs_qp , cost , exitflag]   = quadprog(hh , ff , a_ineq , b_ineq , a_eq , b_eq , [] , [] , [] , options);
if exitflag <= 0
    error('The sign-constrained quadratic fit failed (solver status %d).', exitflag);
end


% finding the binary image of the current AC estimate
binary_image_qp = ac_image(poly_coefs_qp , x_range , y_range);

% computing the associated pixels
pixelized_qp    = pixelization(binary_image_qp , x_range , y_range , spline_degree);

if 1
    p_snr       = 10 * log10( mean(mean(pixelized_noisy.^2)) / mean(mean( (pixelized_qp-pixelized_noisy).^2 )));
    disp(['      >>> pixel consistency: ' , num2str(p_snr) , ' dB , maximum error = ' , num2str(max(max( abs(pixelized_qp-pixelized_noisy) )))])
end






%%%%%%%%%%%%%%
% Checking the qp result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting the difference of the output with the original image
figure,
imshow(abs(binary_image_qp - binary_image))
title('Difference between the reconstructed (zero-Thresh QP) and original images')
xlabel(['PSNR = ' , num2str( -10*log10( mean(mean((binary_image_qp - binary_image).^2)) ) )])

drawnow
