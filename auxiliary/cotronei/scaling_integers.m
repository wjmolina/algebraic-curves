function [f]=scaling_integers(c)
%Computes  the scaling function at integers from the mask c

n_total=length(c)-1; %supporto=[0,n_total]
c=2*c/sum(c); %normalize to sum 2
for i=2:n_total
    for j=2:n_total
        if 2*i-j>=1 & 2*i-j<=n_total+1
            a_matrix(i-1,j-1)=c(2*i-j);
        end
    end
end
[v_matrix,d_array] = eig(a_matrix);

[d,ind] = sort(diag(d_array));
ds = d_array(ind,ind);
vs = v_matrix(:,ind);

f=vs(:,end);
f=[0;f;0]/sum(f); %valori normalizzati a 1 (partition of unity)

    


