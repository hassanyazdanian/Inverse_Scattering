function result = ite_sweep(dimension,N,candidates,reference_k,quadrature_order)
%ITE_SWEEP Match the first computed wavenumber over a candidate grid.
if nargin < 5, quadrature_order = 160; end
validateattributes(candidates, {'numeric'}, {'vector','real','finite','positive','nonempty'});
validateattributes(reference_k, {'numeric'}, {'scalar','real','finite','positive'});
matrices = ite_galerkin_matrices(dimension,N,quadrature_order);
first_k = zeros(size(candidates)); residual = first_k;
for j = 1:numel(candidates)
    [k,r] = ite_wavenumbers(matrices,candidates(j));
    first_k(j) = k(1); residual(j) = r(1);
end
error_abs = abs(first_k-reference_k); [~,best] = min(error_abs);
result = struct('n',candidates,'first_k',first_k,'reference_k',reference_k, ...
    'absolute_error',error_abs,'relative_error',error_abs/reference_k, ...
    'estimated_n',candidates(best),'residual',residual,'dimension',dimension, ...
    'N',N,'quadrature_order',quadrature_order);
end
