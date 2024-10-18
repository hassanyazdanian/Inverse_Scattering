function A = A1_function_3D (del2_phi,del2_phi_c, n, N)

syms r theta
% rr = [0,1];
% thetaa=[0, 2*pi]; 


for i=1:N
    for j=i
        f(r,theta) = del2_phi(i)*del2_phi_c(j)*r^2*sin(theta);
        fun = matlabFunction(simplify(f,'Steps',4));
%         fun = @(theta,r) (1/(n-1))*(matlabFunction(del2_phi(i).*del2_phi_c(j))).*r;
%         fun1 = @(theta,r) (1/(n-1))*(eval(del2_phi(i).*del2_phi_c(j))).*r;
%         A(i,j) = integral2(fun1,0,2*pi,0,1);  

%         fun =matlabFunction(del2_phi(i)*del2_phi_c(j)*r);
        
        A(i,j) = 2*pi*(1/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated'); %2pi for \phi in spherical coordinate 
%        %exclude origin r=0; 
        % A(i,j) = quad2d(fun,0,2*pi,0,1,'AbsTol', 1e-9, 'RelTol', 1e-9);
    end
end