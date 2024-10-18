function B = B1_function_3D (phi,del2_phi,del2_phi_c, n, N)


syms r theta

B = zeros(N,N);

for i=1:N
    for j=i
        f(r,theta) = del2_phi(i)*phi(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = 2*pi*(n/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');
        %exclude origin r=0; 
        f(r,theta) = phi(i)*del2_phi_c(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = B(i,j) + 2*pi*(1/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');
        %exclude origin r=0; 
    end
end

for i=1:N
    if i==1 
        for j=4
            f(r,theta) = del2_phi(i)*phi(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = 2*pi*(n/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');

        f(r,theta) = phi(i)*del2_phi_c(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = B(i,j) + 2*pi*(1/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');
            B(j,i)= B(i,j);
        end
    end
    
    if i==2
        for j=6
           f(r,theta) = del2_phi(i)*phi(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = 2*pi*(n/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');

        f(r,theta) = phi(i)*del2_phi_c(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = B(i,j) + 2*pi*(1/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');
            B(j,i)= B(i,j);
        end
    end

    if i==3
       for j=8
            f(r,theta) = del2_phi(i)*phi(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = 2*pi*(n/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');

        f(r,theta) = phi(i)*del2_phi_c(j)*r^2*sin(theta);
        fun = matlabFunction(f);
        B(i,j) = B(i,j) + 2*pi*(1/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');
            B(j,i)= B(i,j);
        end 
    end

  
end
        
