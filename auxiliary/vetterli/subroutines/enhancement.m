




%%%%%%%%%%%%%%%%%%
% Iterations to improve consistency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_iterations         = 20;
delta           = 1e-4;     % dx for approximating the gradient

current_err     = 1;
new_err         = 0;
iter            = 0;

% main iterations
while ( ((current_err - new_err) / current_err) > 0.001 ) && ( iter < max_iterations )

    iter        = iter+1;
    
    % the gradient of the coefficient with respect to the pixels
    gradient_matrix         = zeros(length(m_x)*length(m_y)  ,  (sum_degree+1) * (sum_degree+2) / 2);
    for coef_ind = 1 : (sum_degree+1) * (sum_degree+2) / 2
        %coef_ind
        
        % marginally changing the "coef_ind"th coefficient
        new_coef            = poly_coefs_qp;
        new_coef(coef_ind)  = new_coef(coef_ind) + delta;
    
        % finding the binary image for the modified coefficients
        binary_image_new    = ac_image(new_coef , x_range , y_range);

        % computing the associated pixels
        pixelized_new       = pixelization(binary_image_new , x_range , y_range , spline_degree);
    

        aux                 = (pixelized_new - pixelized_qp) / delta;
        gradient_matrix(: , coef_ind)   = aux(:);
        
    end
        
    
    % evaluating the pixelization consistency
    diff    = pixelized_qp - pixelized_noisy;
    diff    = diff(:);
    
    current_err     = norm(diff);
   
    
    % possible modification in coefficients
    change  = - pinv(gradient_matrix) * diff;
    
    
    
    
    %%%%%%%%%%%%%
    % tuning the step size
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % finding the binary image of the new AC estimate
    binary_image_qp = ac_image(poly_coefs_qp + change , x_range , y_range);
    
    % computing the associated pixels
    pixelized_qp    = pixelization(binary_image_qp , x_range , y_range , spline_degree);
    
    % evaluating the new error
    new_err = sum(sum( (pixelized_qp-pixelized_noisy).^2 )).^0.5;
    
    % decreasing the step size if the error has increased
    while new_err > current_err * 1.001
        change  = change / 2;
        
        % finding the binary image of the new AC estimate
        binary_image_qp = ac_image(poly_coefs_qp + change , x_range , y_range);
        
        % computing the associated pixels
        pixelized_qp    = pixelization(binary_image_qp , x_range , y_range , spline_degree);
        
        % evaluating the new error
        new_err = sum(sum( (pixelized_qp-pixelized_noisy).^2 )).^0.5;
    end
    
    % updating the coefficients
    poly_coefs_qp   = poly_coefs_qp + change;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%
    
    
    
    if 1
        p_snr       = 10 * log10( mean(mean(pixelized_noisy.^2)) / mean(mean( (pixelized_qp-pixelized_noisy).^2 )));
        disp(['      >>> pixel consistency: ' , num2str(p_snr) , ' dB , maximum error = ' , num2str(max(max( abs(pixelized_qp-pixelized_noisy) )))])
    end
    
    
end









%%%%%%%%%%%%%%
% Checking the final result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting the difference of the output with the original image
figure,
imshow(abs(binary_image_qp - binary_image))
title('Difference between the reconstructed (QP+gradient) and original images')
xlabel(['PSNR = ' , num2str( -10*log10( mean(mean((binary_image_qp - binary_image).^2)) ) )])





%%%%%%%%%%%%
% Saving the final result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imwrite(binary_image_qp , fullfile(results_dir, ['enhanced_qp_binary_image_snr_' , lower(strrep(num2str(pixel_snr), '.', '_')) , '.tif']))
