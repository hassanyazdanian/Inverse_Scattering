function value = spherical_jn(m,x)
% Spherical Bessel function with its finite value at the origin.
validateattributes(m,{'numeric'},{'scalar','integer','nonnegative'});
if isa(x,'sym')
    value = sqrt(pi./(2*x)).*besselj(m+sym(1)/2,x); return
end
value = zeros(size(x)); nonzero = x~=0;
value(nonzero) = sqrt(pi./(2*x(nonzero))).*besselj(m+0.5,x(nonzero));
value(~nonzero) = (m==0);
end
