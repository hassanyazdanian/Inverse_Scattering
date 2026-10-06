% This bracket selects the m=0 root; it is not a global search over orders.
addpath(fileparts(mfilename('fullpath')));
n = 40; m = 0;
k = fzero(@(k) f(k,n,m),[0.70 0.73]);
fprintf('n=%g, m=%d, k=%.15g\n',n,m,k);
