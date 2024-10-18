function y = f_basis_eig(k,r,m)



% J = besselj(m,k*r); % Bessel function of first kind
% dJ = (m*besselj(m, k*r))/r - k*besselj(m + 1, k*r);
% 
% I  = besseli(m, k*r);
% dI = k*besseli(m + 1, k*r) + (m*besseli(m, k*r))/r;

% y = J*dI - dJ*I;
y = k*(besselj(m,k*r)*besseli(m + 1, k*r)+besseli(m, k*r)*besselj(m + 1, k*r));