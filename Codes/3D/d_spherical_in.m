function din = d_spherical_in(nu, x)
% returns the derivative of modifiedspherical Bessel functions 
% of the first kind jnu(x)
% x and nu are scalars

din = spherical_in(nu-1, x)-(nu+1)/x*spherical_in(nu, x);