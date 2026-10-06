% Constant-index unit disk reconstruction. Run from any directory.
script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir,'..','common'));
N = 20; candidates = 2:0.1:20.3; reference_k = 4.159235819877813;
result = ite_sweep(2,N,candidates,reference_k);
fprintf('Estimated n = %.8g (N=%d)\n',result.estimated_n,N);
figure; plot(result.n,result.relative_error,'LineWidth',1.5); grid on;
xlabel('coefficient n'); ylabel('|k_0^N-k_0|/k_0');
title('2D constant-index reconstruction');
output_dir = fullfile(script_dir,'..','..','outputs');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
save(fullfile(output_dir,'reconstruction_2d.mat'),'result','-v7');
print(fullfile(output_dir,'reconstruction_2d.png'),'-dpng','-r150');
