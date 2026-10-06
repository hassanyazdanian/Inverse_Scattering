function matrix = ite_symbolic_matrix(kind,dimension,phi,lap,n,N)
%ITE_SYMBOLIC_MATRIX Optional slow reference assembly (Symbolic Toolbox).
validateattributes(n, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(N, {'numeric'}, {'scalar','integer','>=',1,'<=',20});
if n == 1, error('n must differ from 1.'); end
syms r theta real
if dimension == 2
    weight = r; angular_upper = 2*pi;
elseif dimension == 3
    weight = 2*pi*r^2*sin(theta); angular_upper = pi;
else
    error('dimension must be 2 or 3.');
end
matrix = zeros(N,N);
for i = 1:N
    for j = 1:N
        switch upper(kind)
            case 'A', expression = lap(i)*conj(lap(j))/(n-1);
            case 'B', expression = (n*lap(i)*conj(phi(j))+phi(i)*conj(lap(j)))/(n-1);
            case 'C', expression = n*phi(i)*conj(phi(j))/(n-1);
            otherwise, error('kind must be A, B, or C.');
        end
        fun = matlabFunction(expression*weight,'Vars',{r,theta});
        matrix(i,j) = integral2(fun,0,1,0,angular_upper,'AbsTol',1e-9,'RelTol',1e-7);
    end
end
end
