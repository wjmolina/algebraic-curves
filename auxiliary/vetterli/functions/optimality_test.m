clc
clear all
close all





% main parameters
spl_deg         = 8;

x_max           = 16.5;
x_min           = -16.5;

max_power       = 8;

lambda          = ones(1,max_power);



% dependent parameters
dx              = 1e-2;
x_range         = x_min : dx : x_max;
        
kernel          = [x_range ; centered_spline(spl_deg , x_range)];



% cost handle
c_dg            = @(coefs) cost_dg(1 , coefs , kernel , lambda , max_power , x_min , x_max);
c_g             = @(coefs) cost_g(coefs , kernel , lambda , max_power , x_min , x_max);





% loading a data set
%load('coefs_n_6_m_6_maxPower_6_xMax_10.5.mat')
load('coefs_n_8_m_10_max_power_8_x_max_16_5.mat')

coef_g_l     = coef_g(2:end , :) / coef_g(2 , coef_g(1 , :) == 0);
coef_g_l     = [coef_g(1 , :)  ;  coef_g_l];

coef_dg_l    = coef_dg(2:end , :) / centered_spline(10 , 0);
coef_dg_l    = [coef_dg(1 , :) ; coef_dg_l];





%%

% optimization for 'g'
[coef_g_o , error_g_o]    = optimal_coefs_g(kernel , lambda , max_power , x_min , x_max);
coef_g_o_modif   = coef_g_o(2:end , :) / coef_g_o(2 , coef_g_o(1 , :) == 0);
coef_g_o_modif   = [coef_g_o(1,:) ; coef_g_o_modif];

% cost for this choice
c_g(coef_g_o_modif)

% cost for the loaded data
c_g(coef_g_l)




%%

% optimization for 'dg'
[coef_dg_o , error_dg_o]  = optimal_coefs_dg(kernel , lambda , max_power , x_min , x_max);

% cost for this choice
c_dg(coef_dg_o)

% cost for the loaded data
c_dg(coef_dg_l)
























