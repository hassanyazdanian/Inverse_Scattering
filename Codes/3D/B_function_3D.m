function matrix = B_function_3D(phi,del2_phi,del2_phi_c,n,N)
% Compatibility wrapper for full symbolic reference assembly.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','common'));
matrix = ite_symbolic_matrix('B',3,phi,del2_phi,n,N);
end
