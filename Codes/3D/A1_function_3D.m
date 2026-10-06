function matrix = A1_function_3D(del2_phi,del2_phi_c,n,N)
% Compatibility wrapper for full symbolic reference assembly.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','common'));
matrix = ite_symbolic_matrix('A',3,[],del2_phi,n,N);
end
