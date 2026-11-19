% This file is to read the right set of moment generating coefficients

if coef_setting_method == 1     % the coefficients are set based on 
                                % "spline_degree" and "sum_degree" automatically 
    
    maximum_power    = ceil(sum_degree * 1.5);
    
    switch spline_degree
        case 2
            load(['coefs_n_2_m_6_max_power_' , num2str(maximum_power) , '_x_max_14.mat']);
            
        case 4
            load(['coefs_n_4_m_7_max_power_' , num2str(maximum_power) , '_x_max_16.mat']);
            
        case 6
            load(['coefs_n_6_m_10_max_power_' , num2str(maximum_power) , '_x_max_16_5.mat']);
            
        case 8
            load(['coefs_n_8_m_10_max_power_' , num2str(maximum_power) , '_x_max_16_5.mat']);
    end
    
    
    
else        % user-defined coefficients
    
    
    filename    = 'coefs_n_8_m_8_max_power_8_x_max_13_5_shft_n_2.mat';
    load(filename);
    
    
end