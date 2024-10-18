clc
syms k r theta m c a1

% simplifying the Laplacian of the basis functions in the spherical coordinate
% it is for checking only
% the basis are not a function of \phi.

% fun1 = spherical_jn(m,k*r)*cos(m*theta)
% fun1 = sqrt(pi /(2*k*r))*besseli(m + 0.5, k*r)*cos(m*theta)
% simplify(diff(diff(fun1,r),r)+2/r*diff(fun1,r)+1/(r^2*sin(theta))*diff(sin(theta)*diff(fun1,theta),theta))

fun2 = spherical_jn(m,k*r)*legendreP(m, cos(theta));
simplify(diff(diff(fun2,r),r)+2/r*diff(fun2,r)+1/(r^2*sin(theta))*diff(sin(theta)*diff(fun2,theta),theta))

% diff(fun1,r)=(spherical_in(m-1, k*r)-(m+1)/(k*r)*spherical_in(m, k*r))*cos(m*theta)
% 
% din = spherical_in(m-1, k*r)-(m+1)/(k*r)*spherical_in(m, k*r);
% 
% 
% (spherical_in(m-2, k*r)-(2*m+1)/(k*r)*spherical_in(m-1, k*r)+((m+1)/(k*r))^2*spherical_in(m, k*r))*cos(m*theta)
% 
% 2/r*diff(fun1,r)= (2/r*spherical_in(m-1, k*r)-2/r*(m+1)/(k*r)*spherical_in(m, k*r))*cos(m*theta)
% 
% diff(fun1,theta) = -m*spherical_in(m, k*r)*sin(m*theta)
% 
% sin(theta)*diff(fun1,theta) = -m*spherical_in(m, k*r)*sin(m*theta)*sin(theta)
% 
% diff(sin(theta)*diff(fun1,theta),theta)) = 
% 
% -m*spherical_in(m, k*r)*(sin(m*theta)*cos(theta)+m*cos(m*theta)*sin(theta))
% 
% 1/(r^2*sin(theta))*diff(sin(theta)*diff(fun1,theta),theta)=
% -m/r^2*spherical_in(m, k*r)*(sin(m*theta)*cot(theta)+m*cos(m*theta))




% spherical_in(m-2, k*r)*cos(m*theta)
% 
% spherical_in(m-1, k*r)*((2*(k-m)-1)/(k*r))*cos(m*theta)
% 
% spherical_in(m, k*r)*((((m+1)^2-(k*m+1)^2-2*k+1)/(k*r)^2)*cos(m*theta)-m/r^2*sin(m*theta)*cot(theta))


% -(2^(1/2)*pi^(1/2)*besseli(m + 1/2, k*r)*(2*m*sin(theta*(m - 1)) + k^2*r^2*sin(theta*(m - 1)) - k^2*r^2*sin(theta*(m + 1))))/(4*k^2*r^4*sin(theta)*(1/(k*r))^(3/2))
% 
% fun2= sqrt(pi /(2* k*r)) * besseli(m + 0.5, k*r)*cos(m*theta);
% simplify(diff(diff(fun2,r),r)+2/r*diff(fun2,r)+1/(r^2*sin(theta))*diff(sin(theta)*diff(fun2,theta),theta))
% 
% % fun2=sp_besselj(m,k*r)*cos(m*theta)
% % simplify(diff(diff(fun2,r),r)+2/r*diff(fun2,r)+1/(r^2*sin(theta))*diff(sin(theta)*diff(fun2,theta),theta))
% % 
% % -(2^(1/2)*pi^(1/2)*besselj(m + 1/2, k*r)*(2*m*sin(theta*(m - 1)) - k^2*r^2*sin(theta*(m - 1)) + k^2*r^2*sin(theta*(m + 1))))/(4*k^2*r^4*sin(theta)*(1/(k*r))^(3/2))
%  
% % fun2=(a1*sp_besselj(m, k*r) + sp_besseli(m, k*r))*cos(m*theta);
% % simplify(diff(diff(fun2,r),r)+2/r*diff(fun2,r)+1/(r^2*sin(theta))*diff(sin(theta)*diff(fun2,theta),theta))
