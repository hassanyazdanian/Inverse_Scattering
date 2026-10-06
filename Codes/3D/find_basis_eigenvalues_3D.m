addpath(fileparts(mfilename('fullpath')));
m = 0; root = fzero(@(k) f_basis_eig_3D(k,m),[3.8 4.1]);
fprintf('First clamped m=0 basis root: %.15g\n',root);
