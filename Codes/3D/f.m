function value = f(k,n,m,radius)
%F Ball transmission determinant for spherical harmonic order m.
if nargin < 4, radius = 1; end
validateattributes(n,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(radius,{'numeric'},{'scalar','real','finite','positive'});
x = k*radius; y = sqrt(n)*x;
value = k.*(sqrt(n)*spherical_jn(m,x).*d_spherical_jn(m,y) ...
    -spherical_jn(m,y).*d_spherical_jn(m,x));
end
