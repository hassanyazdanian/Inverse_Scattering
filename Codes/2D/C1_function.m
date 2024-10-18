function C = C1_function (phi, n, N)


syms r theta

C = zeros(N,N);
for i=1:N
    for j=i
        f(r,theta) = phi(i)*phi(j)*r;
        fun = matlabFunction(f);
        C(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    end
end
