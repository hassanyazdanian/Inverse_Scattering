function [k,residuals] = ite_wavenumbers(matrices,n)
%ITE_WAVENUMBERS Positive near-real roots of A + k^2 B + k^4 C.
validateattributes(n, {'numeric'}, {'scalar','real','finite','positive'});
if n == 1, error('n=1 has no contrast; the transmission problem degenerates.'); end
% Multiplication by (n-1) removes the common denominator in the weak form.
A = matrices.S; B = n*matrices.G+matrices.G'; C = n*matrices.M;
[vectors,values] = polyeig(A,B,C);
keep = isfinite(values) & real(values)>0 & ...
    abs(imag(values)) <= 1e-8*max(1,abs(values));
values = real(values(keep)); vectors = vectors(:,keep);
[values,idx] = sort(values); vectors = vectors(:,idx);
if isempty(values), error('No positive near-real eigenvalues for this basis.'); end
k = sqrt(values); residuals = zeros(size(k));
for j = 1:numel(k)
    lambda = values(j); v = vectors(:,j);
    residuals(j) = norm((A+lambda*B+lambda^2*C)*v) / ...
        ((norm(A)+lambda*norm(B)+lambda^2*norm(C))*norm(v));
end
end
