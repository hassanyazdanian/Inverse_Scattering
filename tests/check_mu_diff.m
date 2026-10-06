function check_mu_diff(toolbox_path)
%CHECK_MU_DIFF Optional comparison of the analytical disk and mu-diff solvers.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Codes','2D','LSM'));
addpath(genpath(toolbox_path));
angles = (0:16)'*2*pi/17;
exact = circular_far_field(2,3,1,angles);
numerical = penetrable_far_field(2,3,[0;0],1,angles,'mu-diff');
relative_error = norm(exact-numerical,'fro')/norm(exact,'fro');
assert(relative_error < 1e-8);
fprintf('PASS mu-diff/analytical disk comparison: relative error %.3e\n',relative_error);
end
