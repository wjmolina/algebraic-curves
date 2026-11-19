function cost = cost_g(coef_g , kernel , lambda , max_power , x_min , x_max)
%
% cost = cost_g(coef_g , kernel , lambda , max_power , x_min , x_max)
% 
% This function evaluates the cost incurred by representing x^i g(x) by 
% linear combinations of shifted kernels. 'g(.) is an unknown function. In 
% ideal case, we would like to have 
%
%               x^i g(x)   ~   \sum_{m} b(i,m) kernel(x-m)
%
% where ~ shows almost equality and 'b(i,m)' are the coeeficients stored in
% "coef_g". As we do not know 'g(.)', we only compare these equations with
% respect to each other:
%
%                                   cost = 
% sum_{i} ||\sum_{m} b(i,m) * k_array(x-m) - x * \sum_{m} b(i-1,m) * k_array(x-m)||_2^2
%
%
%
%
% "coef_g" is a 2D matrix containing the coefficients b(i,m). This matrix 
%          includes "max_power+1" number of rows, the first of which
%          corresponds to values of 'm' and the following rows to monomial
%          degrees i = 1 ... max_power.
%
%
% "kernel" is a 2*l matrix that represents the kernel function. Its
%          strcuture is given by
%                        _                                           _ 
%                       |   x_1        x_2      x_3     ...    x_l    |
%              kernel = |                                             |
%                       | Ker(x_1)  Ker(x_2)  Ker(x_3)  ...  Ker(x_l) |
%                       *-                                           -*
%          We assume that {x_i}s form an arithmetic progression. It is
%          further assumed that the kernel vanishes beyond the given range.
% 
%
% "lambda" is a weight vector with length "max_power".
%
%
% "max_power" is the maximum degree 'i' of the monomials multiplied by 
%             'g(x)'. Also, the minimum of 'i' is 1 and not 0. 
%
% "x_min" and "x_max" are the bounds on the interval for which we woul like
%         our almost equality to hold. Essentially, we do not care about 
%         the approximation quality outside of this range. These two values
%         should be integers.                     


% primary parameters
dx              = kernel(1,2) - kernel(1,1);            % sampling resolution
dx              = round(dx * 1e7) / 1e7;

x_range         = x_min : dx : x_max;                   % considered range of 'x' 

m_min           = coef_g(1,1);                          % minimum index 'm' for b(i,m)
m_max           = coef_g(1,end);                        % maximum index 'm' for b(i,m)
m_range         = coef_g(1,:);                          % all values of 'm' for b(i,m)
mmm             = length(m_range);




% storing the shifted extended kernels k_array(x-m)
x_extended      = x_min + m_min : dx : x_max + m_max;
x_large         = x_min + 2*m_min : dx : x_max + 2*m_max;

kernel_large    = zeros(size(x_large));
kernel_large( ( x_large >= kernel(1,1)-1e-10 )&( x_large <= kernel(1,end)+1e-10 ) )     = kernel(2,:);

ker_shift       = zeros(mmm  ,  length(x_extended));
for m_ind = 1 : mmm
    aux         = (x_large + m_range(m_ind) >= x_min+m_min-1e-10)&(x_large + m_range(m_ind) <= x_max+m_max+1e-10);
    ker_shift(m_ind , :)    = kernel_large(aux);
end

aux             = (x_extended >= x_min-1e-10)&(x_extended <= x_max+1e-10);
ker_shift       = ker_shift(: , aux);





% implementing the cost function
cost            = lambda(1) * sum( ((coef_g(2 , :) * ker_shift) .* x_range - coef_g(3 , :) * ker_shift).^2 ) * dx;
for pp = 2 : max_power
    cost        = cost + ...
                  lambda(pp) * sum( ( (coef_g(pp+1 , :) * ker_shift) .* x_range - coef_g(pp+2 , :) * ker_shift).^2 ) * dx;
end

