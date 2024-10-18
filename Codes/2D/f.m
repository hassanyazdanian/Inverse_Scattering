function y = f(k,r,n,m)



J1 = besselj(m,k*r); % Bessel function of first kind
dJ1 = (m*besselj(m, k*r))/r - k*besselj(m + 1, k*r);

J2  = besselj(m, k*n^(1/2)*r);
dJ2 = (m*besselj(m, k*n^(1/2)*r))/r - k*n^(1/2)*besselj(m + 1, k*n^(1/2)*r);

y = J1* dJ2- J2*dJ1;