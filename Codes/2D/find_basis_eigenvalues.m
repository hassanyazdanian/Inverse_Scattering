addpath(fileparts(mfilename('fullpath')));
m = 0; root = fzero(@(k) f_basis_eig(k,m),[3 3.4]);
fprintf('First clamped m=0 basis root: %.15g\n',root);
