function djn = d_spherical_jn(nu, x)
% returns the derivative of spherical Bessel functions of the 
% first kind jnu(x)
% x and nu are scalars

% djn = spherical_jn(nu-1, x)-(nu+1)/x*spherical_jn(nu, x);
djn = nu/x*spherical_jn(nu, x)-spherical_jn(nu+1, x);