function value = d_spherical_jn(m,x)
% Derivative with respect to the complete argument x.
validateattributes(m,{'numeric'},{'scalar','integer','nonnegative'});
if isa(x,'sym')
    value = m./x.*spherical_jn(m,x)-spherical_jn(m+1,x); return
end
value = zeros(size(x)); nonzero = x~=0;
value(nonzero) = m./x(nonzero).*spherical_jn(m,x(nonzero))-spherical_jn(m+1,x(nonzero));
value(~nonzero) = (m==1)/3;
end
