function F = penetrable_far_field(k,n,centers,radii,angles,backend)
%PENETRABLE_FAR_FIELD Single analytic disk or multiple disks via mu-diff.
if nargin < 6, backend = 'analytic'; end
if size(centers,1) ~= 2 || numel(radii) ~= size(centers,2)
    error('centers must be 2-by-M, with one radius per disk.');
end
if strcmp(backend,'analytic')
    if numel(radii) ~= 1, error('analytic backend supports one disk.'); end
    F = circular_far_field(k,n,radii(1),angles,false,centers(:,1));
elseif strcmp(backend,'mu-diff')
    if exist('PenetrableSolver','file') ~= 2 || exist('PenetrableFarField','file') ~= 2
        error('Install mu-diff and add its subdirectories with addpath(genpath(...)).');
    end
    angles = angles(:)'; F = complex(zeros(numel(angles)));
    interior_k = k*sqrt(n)*ones(numel(radii),1);
    for j = 1:numel(angles)
        solution = PenetrableSolver(centers,radii,k,interior_k,'PlaneWave',angles(j),'Direct');
        values = PenetrableFarField(solution,angles);
        F(:,j) = values(:);
    end
else
    error('backend must be analytic or mu-diff.');
end
end
