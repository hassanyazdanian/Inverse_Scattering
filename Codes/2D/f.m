function value = f(k,n,m,radius)
%F Disk transmission determinant for angular order m.
if nargin < 4, radius = 1; end
validateattributes(n,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(m,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(radius,{'numeric'},{'scalar','real','finite','positive'});
x = k*radius; y = sqrt(n)*x;
dJx = (besselj(m-1,x)-besselj(m+1,x))/2;
dJy = (besselj(m-1,y)-besselj(m+1,y))/2;
value = k.*(sqrt(n)*besselj(m,x).*dJy-besselj(m,y).*dJx);
end
