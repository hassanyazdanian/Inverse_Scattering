function run_matlab_tests()
%RUN_MATLAB_TESTS Numerical checks runnable with MATLAB or GNU Octave.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'Codes','common'));
disk_path = fullfile(root,'Codes','2D'); ball_path = fullfile(root,'Codes','3D');
addpath(disk_path); addpath(fullfile(disk_path,'LSM'));
[x,w] = gauss_legendre(12,0,1);
for power = 0:20, assert(abs(sum(w.*x.^power)-1/(power+1)) < 2e-14); end
fprintf('PASS Gauss-Legendre moments\n');
assert(abs(fzero(@(k) f(k,3,0),[4.1 4.2])-4.159235819877813) < 1e-11);
assert(abs(fzero(@(k) f_basis_eig(k,0),[3 3.4])-3.196220616582541) < 1e-11);
fprintf('PASS disk analytical roots\n');
rmpath(disk_path); addpath(ball_path); clear f;
assert(abs(fzero(@(k) f(k,40,0),[0.70 0.73])-0.716709033375646) < 1e-11);
assert(abs(fzero(@(k) f_basis_eig_3D(k,0),[3.8 4.1])-3.926602312047923) < 1e-11);
assert(spherical_jn(0,0)==1 && spherical_jn(2,0)==0);
assert(spherical_in(0,0)==1 && spherical_in(2,0)==0);
assert(d_spherical_jn(1,0)==1/3 && d_spherical_in(1,0)==1/3);
assert(d_spherical_jn(0,0)==0 && d_spherical_in(0,0)==0);
fprintf('PASS ball roots and origin limits\n');
for dimension = [2 3]
    [orders,roots] = ite_basis_data(dimension,20);
    for mode = 1:20
        m = orders(mode); z = roots(mode);
        if dimension == 2
            J = besselj(m,z); I = besseli(m,z);
            Jp = (besselj(m-1,z)-besselj(m+1,z))/2;
            Ip = (besseli(m-1,z)+besseli(m+1,z))/2;
        else
            J = spherical_jn(m,z); I = spherical_in(m,z);
            Jp = d_spherical_jn(m,z); Ip = d_spherical_in(m,z);
        end
        % R(1)=0 by construction; verify the clamped derivative too.
        assert(abs(Jp-(J/I)*Ip)/max([abs(Jp),abs((J/I)*Ip),eps]) < 1e-10);
    end
    matrices = ite_galerkin_matrices(dimension,20,160);
    coarse = ite_galerkin_matrices(dimension,20,80);
    assert(norm(matrices.M-eye(20),'fro') < 1e-10);
    assert(norm(matrices.G-matrices.G','fro')/norm(matrices.G,'fro') < 1e-11);
    assert(min(eig(matrices.S)) > 0);
    assert(norm(matrices.S-coarse.S,'fro')/norm(matrices.S,'fro') < 1e-10);
    [k,residual] = ite_wavenumbers(matrices,3);
    assert(all(k>0) && all(isfinite(k)) && max(residual) < 1e-10);
end
fprintf('PASS matrix orthogonality, symmetry, quadrature convergence and QEP residuals\n');
disk = ite_sweep(2,20,2:0.1:20.3,4.159235819877813);
ball = ite_sweep(3,8,3:100,0.716709033375646);
assert(abs(disk.estimated_n-3) < 1e-12 && ball.estimated_n==40);
assert(abs(disk.first_k(11)-4.167392441300) < 1e-8);
assert(abs(ball.first_k(38)-0.717844989271) < 1e-8);
exact = 0.716709033375646;
errors = zeros(1,3); sizes = [8 15 20];
for j = 1:3
    k = ite_wavenumbers(ite_galerkin_matrices(3,sizes(j)),40);
    errors(j) = abs(k(1)-exact);
end
assert(all(diff(errors)<0));
fprintf('PASS both reconstruction sweeps and 3D basis convergence\n');
failed = false; try, ite_galerkin_matrices(2,21); catch, failed = true; end; assert(failed);
failed = false; try, ite_wavenumbers(ite_galerkin_matrices(2,4),1); catch, failed = true; end; assert(failed);
fprintf('PASS unsupported basis size and zero-contrast guards\n');
rng(21); A = randn(9,5)+1i*randn(9,5); rhs = randn(9,3)+1i*randn(9,3);
[U,S,V] = svd(A,'econ'); lambda = 0.13;
actual = tikhonov_svd(U,diag(S),V,rhs,lambda);
expected = (A'*A+lambda^2*eye(5))\(A'*rhs);
assert(norm(actual-expected,'fro')/norm(expected,'fro') < 1e-12);
fprintf('PASS complex Tikhonov normal equations\n');
angles = (0:16)'*2*pi/17;
F = circular_far_field(2,1,1,angles); assert(norm(F,'fro') < 1e-12);
F = circular_far_field(2,3,1,angles); assert(norm(F-F.','fro') < 1e-11);
shift = [0.3;-0.2]; projection = [cos(angles) sin(angles)]*shift;
translated = circular_far_field(2,3,1,angles,false,shift);
expected = F.*exp(2i*bsxfun(@minus,projection',projection));
assert(norm(translated-expected,'fro') < 1e-12);
scan = lsm_eigenvalues(10,[1.2 1.3 1.4],'analytic',false);
shape = lsm_shape(4,2,'analytic',[0;0],1,false);
assert(all(isfinite(scan.indicator)) && all(isfinite(shape.indicator(:))));
fprintf('PASS disk scattering and LSM smoke tests\n');
old = load(fullfile(root,'tests','data','archived_2d_sweep.mat'));
change = max(abs(disk.first_k-old.k0_N)./abs(old.k0_N));
fprintf('Historical sweep maximum relative change: %.6f%% (comparison, not accuracy target)\n',100*change);
fprintf('ALL MATLAB/OCTAVE TESTS PASSED\n');
end
