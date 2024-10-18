function sp_del2=sp_laplacian(phi)

% Laplacian (of basis functions) in polar coordinate 
syms r theta %phi
dr2 = diff(diff(phi,r),r);
dr = 2/r*diff(phi,r);
dtheta2 = 1/(r^2*sin(theta))*diff(sin(theta)*diff(phi,theta),theta);

sp_del2 = dr2 + dr + dtheta2;