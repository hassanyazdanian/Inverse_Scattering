function C = C1_function_3D (phi, n, N)

syms r theta
% for i=1:N
%     for j=1:N
%         C(i,j)= n/(n-1)*int(int(r*phi(i)*conj(phi(j)),theta,thetaa),r,rr);
%     end
% end

for i=1:N
    for j=i
%         fun = @(theta,r) (n/(n-1))*(eval(phi(i).*conj(phi(j)))).*r;
% %         fun = @(theta,r) (n/(n-1))*(eval(phi(i).*phi(j))).*r;
%         C(i,j) = integral2(fun,0,2*pi,0,1);
        f(r,theta) = phi(i)*phi(j)*r^2*sin(theta);
        fun = matlabFunction(simplify(f,'Steps',4));
        C(i,j) = 2*pi*(n/(n-1))*integral2(fun,0.02,1,0,pi,'AbsTol',1e-9,'Method','iterated');
        %exclude origin r=0; 
    end
end