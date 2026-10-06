function solution = tikhonov_svd(U,s,V,rhs,regularization)
%TIKHONOV_SVD Complex Tikhonov solve using MATLAB's A=U*diag(s)*V'.
validateattributes(regularization, {'numeric'}, {'scalar','real','finite','positive'});
s = s(:); q = numel(s);
solution = V(:,1:q)*bsxfun(@times,s./(s.^2+regularization^2),U(:,1:q)'*rhs);
end
