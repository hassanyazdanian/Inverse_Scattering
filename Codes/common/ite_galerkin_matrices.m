function matrices = ite_galerkin_matrices(dimension,N,quadrature_order)
%ITE_GALERKIN_MATRICES Normalized clamped basis, full radial integration.
% Only constant coefficients on a unit disk/ball are supported.
if nargin < 3, quadrature_order = 160; end
[orders,roots] = ite_basis_data(dimension,N);
[r,w] = gauss_legendre(quadrature_order,0,1);
weights = w .* r.^(dimension-1);
R = zeros(numel(r),N); L = R; angular_norm = zeros(1,N);
for j = 1:N
    m = orders(j); z = roots(j);
    if dimension == 2
        J = besselj(m,z*r); I = besseli(m,z*r);
        ratio = besselj(m,z)/besseli(m,z);
        angular_norm(j) = pi*(1+(m==0));
    else
        factor = sqrt(pi./(2*z*r));
        J = factor.*besselj(m+0.5,z*r); I = factor.*besseli(m+0.5,z*r);
        ratio = besselj(m+0.5,z)/besseli(m+0.5,z);
        angular_norm(j) = 4*pi/(2*m+1);
    end
    R(:,j) = J-ratio*I; L(:,j) = -z^2*(J+ratio*I);
    scale = sqrt(sum(weights.*R(:,j).^2)*angular_norm(j));
    R(:,j) = R(:,j)/scale; L(:,j) = L(:,j)/scale;
end
angular = bsxfun(@eq,orders',orders).*repmat(angular_norm,N,1);
S = (L'*bsxfun(@times,weights,L)).*angular;
G = (L'*bsxfun(@times,weights,R)).*angular;
M = (R'*bsxfun(@times,weights,R)).*angular;
matrices = struct('S',(S+S')/2,'G',G,'M',(M+M')/2, ...
    'dimension',dimension,'N',N,'quadrature_order',quadrature_order);
end
