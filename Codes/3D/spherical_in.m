function in = spherical_in(nu, x)
% returns the Modified spherical Bessel functions of the first kind ynu(x)
% x is a vector or it may be a matrix if nu is a scalar
% if nu is a row and x a column vector, the output js is a matrix
% [nnu lnu] = size(nu);
% [nx lx] = size(x);
% xm = repmat(x, 1, lnu);
% is = sqrt(pi ./(2* xm)) .* besseli(nu + 0.5, x);
in = sqrt(pi /(2*x))*besseli(nu + 0.5, x);