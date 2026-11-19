






%%%%%%%%%%%%%%
% pixelizing and all-white image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Generating the pixels of an all-white image ----')
tic

% Applying the spline PSF on the binary image, followed by sampling
white_pixs      = pixelization(ones(length(y_range) , length(x_range)) , x_range , y_range , spline_degree);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')




%%%%%%%%%%%%%%%%%%%%%%%%%
% finding the moments of an all-white image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Generalized moments of an all-white image ----')
tic

% Estimating the generalized moments of the image
[aux1 , aux2 , white_moments]   = moment_generator(white_pixs , m_x , m_y , coef_g , coef_dg);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forming the inequalities regarding the polynomial coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Forming the inequality matrix ----')
tic

% Forming the matrix storing the inequality coefficients
even_moment_matrix  = matrix4inequality(sum_degree , moments_x_y , white_moments);

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')








%%%%%%%%%%%
% quadratic programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Solving the quadratic programming ----')
tic

% initializing the minimization with a centered circle
x0 = zeros(size(moment_matrix, 2), 1);
x0(1) = 36;
if sum_degree >= 2
    x0([4, 6]) = -1;
end

% minimization settings
options         = optimset('MaxIter' , 1e8);

% minimization problem with inequality and equality constraints
[poly_coefs_qp, ~, exit_flag] = quadprog(moment_matrix' * moment_matrix  ,  zeros((sum_degree+1)*(sum_degree+2)/2,1)  ,  ...
                           even_moment_matrix  ,  zeros(size(even_moment_matrix,1),1)  ,  ...
                           [1 , zeros(1,(sum_degree+1)*(sum_degree+2)/2-1)]  ,  [1] , ...
                           [] , [] , x0 , options);
if exit_flag <= 0
    error('The moment-constrained quadratic fit failed (solver status %d).', exit_flag);
end                    
                       
ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')




%%%%%%%%%%%%%%%%
% forming the estimated curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Forming the estimated AC image ----')
tic

binary_image_qp = ac_image(poly_coefs_qp , x_range , y_range);

figure,
imshow(binary_image_qp)

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')



%%%%%%%%%%
% Saving the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Saving the estimated AC image ----')
tic

imwrite(binary_image_qp , fullfile(results_dir, ['qp_binary_image_snr_' , lower(strrep(num2str(pixel_snr), '.', '_')) , '.tif']))

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')