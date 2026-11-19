% The purpose of this subroutine is to read the polynlmial coefficients 
% from a file 

row_ind         = 3;

load('bounded_curve.mat')

poly_coefs      = -c_array(row_ind , :);