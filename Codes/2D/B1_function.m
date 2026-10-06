function matrix = B1_function(phi,del2_phi,del2_phi_c,n,N)
% Compatibility wrapper for full symbolic reference assembly.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','common'));
matrix = ite_symbolic_matrix('B',2,phi,del2_phi,n,N);
end
