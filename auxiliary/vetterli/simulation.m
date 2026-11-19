% This file first generates an image based on the roots of a two-variate
% polynomial; then, converts the image into a number of pixels using a
% given PSF. Finally, the goal is to recover the polynomial based on the
% noisy pixels.

clc
clear all
close all


% Adding the path to the location of supporting files
simulation_root = fileparts(mfilename('fullpath'));
results_dir = fullfile(simulation_root, 'results');
addpath(fullfile(simulation_root, 'functions'), fullfile(simulation_root, 'moment_coefs'), ...
        fullfile(simulation_root, 'poly_coefs'), fullfile(simulation_root, 'subroutines'));




%%
%%%%%%%%%%%%%
% Setting the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('      --- Setting the Parameters ----')
tic

x_range         = [-15:0.01:15];        % the range x_grid values for evaluating the polynomial
y_range         = [-15:0.01:15];        % the range y_grid values for evaluating the polynomial

spline_degree   = 4;        % The PSF is assumed to be the kronecker product of two 
                            % identical B-splines of this degree
                            
sum_degree      = 4;        % assumption about the maximum sum degree in the reconstruction

peak_snr            = 30;       % The image will be subject to additive Gaussian noise with
                            % approaximately this peak_snr value in dB.

% The coefficients of the original two-variate polynomial; the ordering is as follows:
% [x^0y^0   x^0y^1   x^1y^0   x^0y^2   x^1y^1   x^2y^0   x^0y^3   x^1y^2 ...]
coef_method     = 1;        % 1: direct input
                            % 2: loading from a file
switch coef_method
    case 1
        poly_coefs      = [13 8 -4 -2 0 -1];
%         poly_coefs      = [10 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0];
        
    case 2
        poly_coef_load
end



ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')






%%
%%%%%%%%%%%%%%%%%%%%%%%%
% Generating the image and its pixelized version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calling the corresponding subroutine
ac_gen_pixelization





%%
%%%%%%%%%%%%%
% Adding Gaussian noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('      --- Adding noise to the AC image ----')
tic

% variance of the Gaussian noise to produce the desired peak_snr
noise_variance      = 10^(-peak_snr / 10);

% noisy image
binary_image_noisy  = imnoise(binary_image , 'gaussian' , 0 , noise_variance);

actual_psnr         = -10 * log10( mean(mean((binary_image - binary_image_noisy).^2)) );
disp(['          The actual PSNR is ' , num2str(actual_psnr) , 'dB'])

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Regenerating noisy pixels ----')
tic

% Applying the spline PSF on the noisy image, followed by sampling
noise_method    = 2;        % 1: pixelizing the noisy binary image
                            % 2: directing adding noise to the pixels
                            
switch noise_method
    case 1
        [pixelized_noisy , m_x , m_y]   = pixelization(binary_image_noisy , x_range , y_range , spline_degree);
        
    case 2
        %noise_im        = randn(size(pixelized));
        %pixelized_noisy = pixelized + 10^(-peak_snr / 20) * sqrt(mean(mean(pixelized.^2)) / mean(mean(noise_im.^2))) * noise_im ;
        pixelized_noisy  = imnoise(pixelized , 'gaussian' , 0 , noise_variance * mean(mean(pixelized.^2)));
        
end

pixel_snr           = 10 * log10( mean(mean(pixelized.^2)) / mean(mean((pixelized - pixelized_noisy).^2)) );
disp(['          The pixel SNR is ' , num2str(pixel_snr) , 'dB'])

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


disp('      --- Saving the results ----')
tic

% Saving the results
imwrite(binary_image_noisy , fullfile(results_dir, ['noisy_binary_image_psnr_' , lower(strrep(num2str(peak_snr), '.', '_')) , '.tif']))
imwrite(imresize(pixelized_noisy , round(length(binary_image_noisy) / length(pixelized_noisy)) , 'nearest') , ...
        fullfile(results_dir, ['noisy_pixelized_snr_' , lower(strrep(num2str(pixel_snr), '.', '_')) , '.tif']))

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')





%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating polynomial coefficients based on noisy pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('      --- Loading the set of coefficients ----')
tic

% choosing the set of coefficient
coef_setting_method     = 1;    % 1: Automatic selection
                                % 2: User-defined (should be set in the file "load_set.m")
                        
                        
% calling the subroutine to read the required file                        
load_set

ttt             = toc;
disp(['          It took ' , num2str(ttt) , ' seconds'])
disp(' ')


%%%%%%%%%%%%%%


% calling the subroutine to perform the rest of the reconstrcution
reconstructor