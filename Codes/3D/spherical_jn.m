function jn = spherical_jn(nu, x)
% returns the spherical Bessel functions of the first kind jnu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix
[nnu lnu] = size(nu);
[nx lx] = size(x);
xm = repmat(x, 1, lnu);
jn = sqrt(pi ./(2* xm)) .* besselj(nu + 0.5, x);
% jn = sqrt(pi/(2*x))*besselj(nu + 0.5, x);