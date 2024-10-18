syms k r theta m c

% simplifying the Laplacian of the basis functions in the polar coordinate
% it is for checking only

fun1=besseli(m,k*r)*cos(m*theta)
simplify(diff(diff(fun1,r),r)+1/r*diff(fun1,r)+1/r^2*diff(diff(fun1,theta),theta))
% k^2*cos(m*theta)*besseli(m, k*r)

fun2=besselj(2,k*r)*cos(2*theta)
simplify(diff(diff(fun2,r),r)+1/r*diff(fun2,r)+1/r^2*diff(diff(fun2,theta),theta))
% -k^2*cos(m*theta)*besselj(m, k*r)

cos(m*theta)*(k*(k*besseli(m, k*r) - (besseli(m + 1, k*r)*(m + 1))/r) - (m*besseli(m, k*r))/r^2 + (m*(k*besseli(m + 1, k*r) + (m*besseli(m, k*r))/r))/r)