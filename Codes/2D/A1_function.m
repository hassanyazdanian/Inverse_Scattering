function A = A1_function (del2_phi,del2_phi_c, n, N)

syms r theta

A = zeros(N,N);
for i=1:N
    for j=i
        f(r,theta) = del2_phi(i)*del2_phi_c(j)*r;
        fun = matlabFunction(f);       
        A(i,j) = (1/(n-1))*integral2(fun,0,1,0,2*pi);
%         A(i,j) = quad2d(fun,0,2*pi,0,1,'AbsTol', 1e-9, 'RelTol', 1e-9);
    end
end
