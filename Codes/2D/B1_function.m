function B = B1_function (phi,del2_phi,del2_phi_c, n, N)

%Warning: This function works ONLY for N = 20!!! 
syms r theta

B = zeros(N,N);

for i=1:N
    for j=i
        f(r,theta) = del2_phi(i)*phi(j)*r;
        fun = matlabFunction(f);
        B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);

        f(r,theta) = phi(i)*del2_phi_c(j)*r;
        fun = matlabFunction(f);
        B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
    end
end

for i=1:N
    if i==1 
        for j=[4 9 17]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end
    
    if i==2
        for j=[6 13]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end

    if i==3
       for j=[8 16]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end 
    end

    if i==4
        for j=[9 17 19]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end

    if i==5
        for j=[11 14 18 20]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end

    if i==6
        for j=13
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end

    if i==7
        for j=[15 18]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end

    if i==8
        for j=16
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 

    if i==9
        for j=[17 19]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 

    if i==10
        for j=[14 17 19]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 

    if i==11
        for j=[14 16 18 20]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 

    if i==12
        for j=14
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 

    if i==14
        for j=[15 18 19]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 

   if i==15
        for j=[18 20]
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
   end 

  if i==16
    for j=[11 20]
        f(r,theta) = del2_phi(i)*phi(j)*r;
        fun = matlabFunction(f);
        B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);

        f(r,theta) = phi(i)*del2_phi_c(j)*r;
        fun = matlabFunction(f);
        B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
        B(j,i)= B(i,j);
    end
  end 

 if i==17
        for j=20
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
 end 
 
   if i==18
        for j=20
            f(r,theta) = del2_phi(i)*phi(j)*r;
            fun = matlabFunction(f);
            B(i,j) = (n/(n-1))*integral2(fun,0,1,0,2*pi);
    
            f(r,theta) = phi(i)*del2_phi_c(j)*r;
            fun = matlabFunction(f);
            B(i,j) = B(i,j) + (1/(n-1))*integral2(fun,0,1,0,2*pi);
            B(j,i)= B(i,j);
        end
    end 
end
        
