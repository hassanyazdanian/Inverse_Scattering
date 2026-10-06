function [x,w] = gauss_legendre(N,a,b)
%GAUSS_LEGENDRE Nodes and weights on [a,b] (Golub-Welsch construction).
validateattributes(N, {'numeric'}, {'scalar','integer','>=',2});
validateattributes(a, {'numeric'}, {'scalar','real','finite'});
validateattributes(b, {'numeric'}, {'scalar','real','finite','>',a});
j = (1:N-1)'; beta = j ./ sqrt(4*j.^2-1);
[V,D] = eig(diag(beta,1)+diag(beta,-1));
[x,idx] = sort(diag(D)); w = 2*(V(1,idx)').^2;
x = (a+b)/2+(b-a)/2*x; w = (b-a)/2*w;
end
