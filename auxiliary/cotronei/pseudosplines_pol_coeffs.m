function [c,offset]=pseudosplines_pol_coeffs(m,k,p,a,b)

% Computes the coefficients associated to the reproduction of monomials of
% degree k (k=0,...,m-1), in the interval [a,b], by the pseudospline with mask p and with reproduction order m
% Note that the filter p is normalized so that the sum of the filter coeffs is 2

n_total=length(p); %filter support [0,n_total-1]. 
offset=(n_total-1)/2; %IMPORTANT: this is the shift tau needed for evaluating correctly the coefficients as values of the monomials
p=2*p/sum(p); % Normalize so that sum=2
m_matrix=n_total-1;
 c=([-(m_matrix-a-1):b-1]+offset).^k;
 