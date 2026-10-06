function result = lsm_eigenvalues(n,k_values,backend,make_plot)
%LSM_EIGENVALUES Regularized spectral indicator for a unit penetrable disk.
% Adapted from the archived project/mu-diff example.
% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry,
% University of Lorraine, CNRS, France. GPL-3.0-or-later.
% See LICENSE-mu-diff.txt and THIRD_PARTY_NOTICES.md.
if nargin < 1, n = 10; end
if nargin < 2, k_values = linspace(0.5,7,261); end
if nargin < 3, backend = 'analytic'; end
if nargin < 4, make_plot = true; end
validateattributes(k_values,{'numeric'},{'vector','positive','finite','real','nonempty'});
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','common'));
count = 61; angles = (0:count-1)'*2*pi/count; weight = 2*pi/count;
regularization = 1e-6; indicator = zeros(size(k_values));
for j = 1:numel(k_values)
    F = penetrable_far_field(k_values(j),n,[0;0],1,angles,backend);
    [U,S,V] = svd(weight*F,'econ');
    solution = tikhonov_svd(U,diag(S),V,sqrt(weight)*ones(count,1),regularization);
    indicator(j) = norm(solution);
end
result = struct('k',k_values,'indicator',indicator,'n',n, ...
    'regularization',regularization,'backend',backend);
if make_plot
    figure; semilogy(k_values,indicator,'LineWidth',1.5); grid on;
    xlabel('Wavenumber k'); ylabel('Regularized solution norm');
    title(sprintf('Unit disk spectral indicator, n=%g',n));
end
end
