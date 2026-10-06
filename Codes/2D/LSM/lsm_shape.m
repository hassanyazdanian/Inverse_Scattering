function result = lsm_shape(n,k,backend,centers,radii,make_plot,sound_soft)
%LSM_SHAPE Linear sampling of disks using an area-weighted far-field matrix.
% Adapted from the archived project/mu-diff example.
% mu-diff - Copyright (C) 2014-2020 X. Antoine and B. Thierry,
% University of Lorraine, CNRS, France. GPL-3.0-or-later.
% See LICENSE-mu-diff.txt and THIRD_PARTY_NOTICES.md.
if nargin < 1, n = 4; end
if nargin < 2, k = 2*pi; end
if nargin < 3, backend = 'analytic'; end
if nargin < 4, centers = [0;0]; end
if nargin < 5, radii = 1; end
if nargin < 6, make_plot = true; end
if nargin < 7, sound_soft = false; end
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','common'));
count = 81; angles = (0:count-1)'*2*pi/count; weight = 2*pi/count;
if sound_soft
    if ~strcmp(backend,'analytic') || numel(radii) ~= 1
        error('sound_soft is available for one analytic disk.');
    end
    F = circular_far_field(k,n,radii,angles,true,centers);
else
    F = penetrable_far_field(k,n,centers,radii,angles,backend);
end
[U,S,V] = svd(weight*F,'econ'); regularization = 1e-5;
extent = max(2,max(abs(centers(:)))+max(radii(:))+0.5);
axis_points = linspace(-extent,extent,101); [X,Y] = meshgrid(axis_points);
points = [X(:)';Y(:)']; directions = [cos(angles) sin(angles)];
indicator = zeros(1,numel(X));
for start = 1:1024:numel(X)
    indices = start:min(start+1023,numel(X));
    rhs = sqrt(weight)*exp(-1i*k*directions*points(:,indices));
    solution = tikhonov_svd(U,diag(S),V,rhs,regularization);
    indicator(indices) = log10(max(sqrt(sum(abs(solution).^2,1)),realmin));
end
indicator = reshape(indicator,size(X));
result = struct('x',axis_points,'y',axis_points,'indicator',indicator, ...
    'n',n,'k',k,'backend',backend,'regularization',regularization);
if make_plot
    figure; imagesc(axis_points,axis_points,indicator); axis xy equal tight; colorbar;
    xlabel('x'); ylabel('y'); title('log10 regularized solution norm'); hold on;
    theta = linspace(0,2*pi,200);
    for j = 1:numel(radii)
        plot(centers(1,j)+radii(j)*cos(theta),centers(2,j)+radii(j)*sin(theta),'w--');
    end
end
end
