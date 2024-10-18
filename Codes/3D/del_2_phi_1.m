function del_2_phi = del_2_phi_1(a,m,k)

syms r theta

% nc= -1.016170554350756e+04; %normalization_coeff

fun = (a*spherical_jn(m,k*r)+spherical_in(m,k*r))*legendreP(m, cos(theta));

% del_2_phi = simplify(diff(diff(fun,r),r)+2/r*diff(fun,r) ...
%     +1/(r^2*sin(theta))*diff(sin(theta)*diff(fun,theta),theta),'Steps',10);
   
del_2_phi = simplify(diff(fun,r,2)+2/r*diff(fun,r) ...
    +1/(r^2*sin(theta))*diff(sin(theta)*diff(fun,theta),theta),'Steps',4);

