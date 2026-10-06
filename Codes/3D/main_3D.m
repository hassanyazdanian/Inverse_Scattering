% Constant-index unit ball reconstruction. Run from any directory.
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'..','common'));
N = 8; candidates = 3:100; reference_k = 0.716709033375646;
result = ite_sweep(3,N,candidates,reference_k);
fprintf('Estimated n = %.8g (N=%d)\n',result.estimated_n,N);
figure; plot(result.n,result.relative_error,'LineWidth',1.5); grid on;
xlabel('coefficient n'); ylabel('|k_0^N-k_0|/k_0');
title('3D constant-index reconstruction');
output_dir = fullfile(script_dir,'..','..','outputs');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
save(fullfile(output_dir,'reconstruction_3d.mat'),'result','-v7');
print(fullfile(output_dir,'reconstruction_3d.png'),'-dpng','-r150');
