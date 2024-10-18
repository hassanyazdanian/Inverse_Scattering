function d_js = d_sp_besselj(nu, x)
% returns the derivative of spherical Bessel functions jnu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix
[nnu lnu] = size(nu);
[nx lx] = size(x);
xm = repmat(x, 1, lnu);
d_js = ((nu)/xm*sphbes(nu, x)- sphbes(nu+1, x));