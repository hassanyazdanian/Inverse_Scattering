function [phi,del2_phi,del2_phi_c] = basis_functions(N)
% Symbolic reference basis; no plotting, exports, or workspace clearing.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','common'));
[orders,roots] = ite_basis_data(2,N);
syms r theta real
phi = sym(zeros(N,1)); del2_phi = phi;
for j = 1:N
    m = orders(j); z = roots(j);

    J = besselj(m,z*r); I = besseli(m,z*r);
    ratio = besselj(m,z)/besseli(m,z); angular = cos(m*theta);

    phi(j) = (J-ratio*I)*angular;
    del2_phi(j) = -z^2*(J+ratio*I)*angular;
end
del2_phi_c = conj(del2_phi);
end
