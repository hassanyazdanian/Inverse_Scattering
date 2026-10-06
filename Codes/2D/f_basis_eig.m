function value = f_basis_eig(k,m)
% Clamped disk characteristic equation, up to an irrelevant sign.
value = k.*(besselj(m,k).*besseli(m+1,k)+besseli(m,k).*besselj(m+1,k));
end
