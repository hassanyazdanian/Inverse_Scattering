function value = f_basis_eig_3D(k,m)
% Clamped ball characteristic equation.
value = spherical_jn(m,k).*d_spherical_in(m,k) ...
    -d_spherical_jn(m,k).*spherical_in(m,k);
end
