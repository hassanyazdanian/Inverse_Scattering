function y = f(k,r,n,m)



j1 = spherical_jn(m,k*r); % Spherical Bessel function of first kind
% dj1 = d_spherical_jn(m,k*r);
dj1 = pi^(1/2)*(k*besselj(m - 1/2, k*r) - (besselj(m + 1/2, k*r)*(m + 1/2))/r)*(1/(2*k*r))^(1/2) - (pi^(1/2)*besselj(m + 1/2, k*r))/(4*k*r^2*(1/(2*k*r))^(1/2));        

j2  = spherical_jn(m, k*n^(1/2)*r); 
% dj2 =  d_spherical_jn(m, k*n^(1/2)*r);
dj2 = - pi^(1/2)*((besselj(m + 1/2, k*n^(1/2)*r)*(m + 1/2))/r - k*n^(1/2)*besselj(m - 1/2, k*n^(1/2)*r))*(1/(2*k*n^(1/2)*r))^(1/2) - (pi^(1/2)*besselj(m + 1/2, k*n^(1/2)*r))/(4*k*n^(1/2)*r^2*(1/(2*k*n^(1/2)*r))^(1/2));       
y = j1* dj2- j2*dj1;