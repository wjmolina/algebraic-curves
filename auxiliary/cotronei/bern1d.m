function b_matrix=bern1d(n,i,x,a,b)
% Computes the i-th Bernstein polynomial of degree n supported in [a,b]

        b_matrix=1/(b-a)^n*nchoosek(n,i)*(x-a).^i.*(b-x).^(n-i);
    
