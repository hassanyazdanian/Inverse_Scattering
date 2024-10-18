function y = f_basis_eig_3D(k,r,m)



j = spherical_jn(m,k*r); % Spherical Bessel function of the first kind
% dj = pi^(1/2)*(k*besselj(m - 1/2, k*r) - (besselj(m + 1/2, k*r)*(m + 1/2))/r)*(1/(2*k*r))^(1/2) - (pi^(1/2)*besselj(m + 1/2, k*r))/(4*k*r^2*(1/(2*k*r))^(1/2));
dj = d_spherical_jn(m,k*r);
% 
i  = spherical_in(m,k*r);% Spherical Bessel function of the second kind
% di = pi^(1/2)*(k*besseli(m - 1/2, k*r) - (besseli(m + 1/2, k*r)*(m + 1/2))/r)*(1/(2*k*r))^(1/2) - (pi^(1/2)*besseli(m + 1/2, k*r))/(4*k*r^2*(1/(2*k*r))^(1/2));
% dI = k*besseli(m + 1, k*r) + (m*besseli(m, k*r))/r;
di = d_spherical_in(m,k*r);

y = j*di - dj*i;
