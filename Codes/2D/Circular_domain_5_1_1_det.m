addpath(fileparts(mfilename('fullpath')));
n = 10; k = linspace(0.01,10,2000);
figure; hold on;
for m = 0:3, plot(k,f(k,n,m),'DisplayName',sprintf('m=%d',m)); end
grid on; xlabel('k'); ylabel('transmission determinant'); legend('show');
title(sprintf('2D determinant, n=%g',n));
